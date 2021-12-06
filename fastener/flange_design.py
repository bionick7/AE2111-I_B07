import numpy as np
import math

"""
Calculates the design of the lug. Also includes mass calculation, since it's used locally as well

4th version
"""

BOUNDARIES = {
    "Minimum thickness": 1,  # As a fraction of Diameter
    "Minimum width": 1.01,  # As a fraction of Diameter
    "Maximum width": 10,  # As a fraction of diameter
    "Number of flanges": 1,
}

ALLOWABLE_BOLT_STRESS = 400e6 / 1.25  # MPa Stainless steel bolt (in shear)


def calculate_flange_dimensions(allowable, forces, previous_design_parameters, t_steps, w_iterations=6, w_steps=20):
    """
    :param allowable: maximum allowable stress
    :param forces: current lug-forces
    :param previous_design_parameters: design parameters from the previous iteration / initial design parameters
    :param t_steps: number of thickness values tested
    :param w_iterations: precision iterations for width
    :param w_steps: number of widths tested per iteration
    :return: thickness, width, mass tuple for the lightest configuration, D1 diameter
    """
    diameter = get_bolt_diameter(forces)  # Set the diameter based on the bolt materials

    backplate_offset = previous_design_parameters["lug x offset"]
    min_width = max(3 * previous_design_parameters["D2"][0], previous_design_parameters["min fastener outer"])

    # Unpack some relevant variables
    t_min = BOUNDARIES["Minimum thickness"] * diameter
    f_z = forces["Fx"]

    t_min = max(math.sqrt(3 * f_z / allowable), 0.05 * diameter, t_min)  # Make sure t_min does not go below a t_over_d < 0.05
    t_max = 25 * t_min
    var_thickness_results = []  # Each item in this list is in the form:
    for thickness in np.linspace(t_min, t_max, t_steps):  # For all thicknesses in the range for which the equations are valid:
        width = get_width(diameter, thickness, allowable, forces, backplate_offset, min_width)
        if width is None:  # If the data point was invalid, ignore it
            continue
        params = previous_design_parameters.copy()
        params["t1"] = thickness
        params["w"] = width
        m = mass(params, 1)  # Density is irrelevant
        var_thickness_results.append((thickness, width, m))  # Append a tuple containing all information about the data point to the list of data points
    var_thickness_results.sort(key=lambda x: x[2])
    best = var_thickness_results[0]
    return best, diameter


def get_width(diameter, thickness, allowable, forces, backplate_offset, min_width, iterations=6, w_steps=20):
    """
    :param diameter: D1 diamter
    :param forces: current forces on the lug
    :param allowable: allowable stress
    :param backplate_offset: additional offset planned for the backplate
    :param min_width: minimum allowable width
    :param iterations: number of iterations. Each iteration increases precision
    :return: The minimum required width for a given thickness
    """
    w_min = max(BOUNDARIES["Minimum width"] * diameter, min_width)
    w_max = BOUNDARIES["Maximum width"] * diameter
    n = BOUNDARIES["Number of flanges"]
    # + the margin

    # Get all forces (and calculate the equivalent horizontal force)
    f_z = forces["Fz"]
    f_x = forces["Fx"]
    f_y = forces["Fy"]
    m_z = forces["Mz"]

    # Calculate the width required to prevent shear out of the pin
    for _ in range(iterations):  # Each iteration reduces the range by a factor 10
        data = []
        step_size = (w_max - w_min) / w_steps
        for width in np.linspace(w_min, w_max, (w_steps + 1)):
            # Get k_ty and k_bry
            k_ty = get_k_ty(width, thickness, diameter)
            k_bry = get_k_bry(width, thickness, diameter)
            if k_ty <= 0 or k_bry <= 0:  # If any constant turns out to be less than 0, the value is outside of the accurate domain.
                continue
            new_diameter = (((f_y * k_ty / n) ** 1.6 + (f_z * k_bry) ** 1.6) / (k_ty * k_bry * allowable * thickness) ** 1.6) ** (1 / 1.6)
            diameter_ratio = abs(math.log(diameter / new_diameter))  # Test how close the value is to the actual
            # (original) value. If it is close, it means the value is close to the convergeance point. Log is used to
            # ensure that both a ratio of 0.5 and 2 receive the same grade.
            data.append([diameter_ratio, width])
        if len(data) == 0:  # If no valid data points were found, return a width of 0 to indicate this.
            return 0
        closest = sorted(data)[0][1]  # Extract the first item from the list of sorted items (the sorted method sorts
        # based on the first item in each sublist, which is the ratio. The lower the ratio, the closer the value.)
        # Both new values for upper and lower bounds are clamped between the same limits as the iteration before,
        # to prevent the values from going outside of the range of the equations for k_\\.
        w_min = max(w_min, closest - step_size)  # Set the new lower bound to the diameter before the closest
        w_max = min(w_max, closest + step_size)  # and the new upper bound to the diameter after the closest.
        # This should always include the closest value, unless the graph makes a very weird shape.
    # If the value is "close to" (= equal to, but taking into account floating point inaccuracy), disregard it,
    # as it most likely will NOT be an accurate value, as this almost certainly means the correct value is outside of the valid region.

    # Calculate the width required to prevent bending in the flange
    if allowable * thickness**2 == 3 * f_z:
        return None

    sideways_bending = (6 * f_x * backplate_offset + 6 * m_z) / abs(allowable * thickness**2 - 3 * f_z)

    a = allowable * thickness
    b = -3 * f_z
    c = -6 * backplate_offset * f_z

    vertical_bending = abs((-b + math.sqrt(b*b-4*a*c))/(2*a))

    #print((6 * (vertical_bending/2 + backplate_offset) * f_y) / (vertical_bending**2 * thickness * sigma_y))
    
    if math.isclose(closest, BOUNDARIES["Minimum width"] * diameter):
        return None
    elif math.isclose(closest, BOUNDARIES["Maximum width"] * diameter):
        return None
    return max(closest, sideways_bending, vertical_bending)


def get_k_ty(width, thickness, diameter):
    """
    Index emperically found k_ty values. parameters found through linear regression
    """
    # Calculate the values for the Area's
    A_av = 6 * thickness / ((8 / (width - diameter)) + (4 / (width - (diameter * math.sqrt(2) / 2))))
    A_br = diameter * thickness
    x = A_av / A_br
    k_ty = 0.2338 * x**4 - 0.6482 * x**3 + 0.2368 * x**2 + 1.2154 * x  # Formula found using excel trendlines
    return k_ty


def get_k_bry(width, thickness, diameter):
    """
    Index emperically found k_bry values. parameters found through linear regression
    """
    t_over_d = thickness / diameter
    t_over_d = min(t_over_d, 0.5)  # If t_over_d > 0.5, use 0.5 as the original source says the line for t_over_d > 0.5 is equal to t_over_d = 0.5
    x = width / (2 * diameter)
    # Calculate the coefficients of the equation (a_n is should be multiplied with x^n).
    a_4 = t_over_d * 0.444 - 0.134
    a_3 = t_over_d * -3.69 + 1.2848
    a_2 = t_over_d * 9.9   - 4.144
    a_1 = t_over_d * -8.06 + 5.807
    a_0 = t_over_d * 1.9   - 2.003
    k_bry = a_4 * x**4 + a_3 * x**3 + a_2 * x**2 + a_1 * x + a_0
    return k_bry


def mass(design_parameters, density, round_corners=False, debug=False):
    """
    :param design_parameters: current design paramter dictionary
    :param density: density of the material
    :param round_corners: if true, rounds the corners of the backplate
    :return: mass of the entire lug
    """
    w = design_parameters["w"]
    surface_lug = (1/2 + np.pi/8) * w * w  # Surface of lug without hole
    surface_lug -= np.pi / 4 * (design_parameters["D1"] ** 2)  # Subtract hole surface
    surface_lug += w * design_parameters["lug x offset"]  # Additional offset
    volume_lug = surface_lug * design_parameters["t1"]  # Multiply whith thickness

    hole_area = sum(design_parameters["D2"] * design_parameters["D2"] / 4 * np.pi)
    base_area = w * design_parameters["x_width"] - hole_area
    if round_corners:
        base_area -= w*w*(1 - np.pi/4)
    if debug:
        print(volume_lug + base_area * design_parameters["t2"])
    return volume_lug * density + base_area * design_parameters["t2"] * density


def get_bolt_diameter(forces):
    """
    :param forces: Current lug-forces
    :return: The minimal diameter of the D1 hole/bolt
    """
    f_z = forces["Fz"]
    f_x = forces["Fx"]
    n = BOUNDARIES["Number of flanges"]
    force = math.sqrt(f_z ** 2 + f_x ** 2) / (2*n)
    area = force / ALLOWABLE_BOLT_STRESS * 4/3
    diameter = 2 * math.sqrt(area / np.pi)
    return diameter
