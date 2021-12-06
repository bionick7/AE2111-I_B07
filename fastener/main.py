import numpy as np

import fastener_forces
import flange_design
import stress_checks
import thermal_stress_check as thermal

"""
    The main program. Will run the design algorithm for all given materials and print out the result.
    Change initial parameters as needed. If FIX_FASTENER_DIAMETER will try to force given FASTENER_DIAMETER,
    unless the forces exceed it. All units in metric
"""

# ================================
# CONSTANT PARAMETERS
# ================================

MATERIALS = {
            # https://www.gabrian.com/2024-aluminum-properties/
            "Aluminium-2024-T4": {
                "Yield strength": 324e6,  # Pa
                "Ultimate strength": 469e6,  # Pa
                "Youngs modulus": 73.1e9,  # Pa
                "Density": 2780,  # kg/m^3
                "Thermal expansion coefficient": 23.2e-6  # 1/K
            },
            # https://www.gabrian.com/7075-aluminum-properties/
            "Aluminium-7075-T6": {
                "Yield strength": 503e6,  # Pa
                "Ultimate strength": 572e6,  # Pa
                "Youngs modulus": 71.7e9,  # Pa
                "Density": 2810,  # kg/m^3
                "Thermal expansion coefficient": 23.4e-6  # 1/K
            },
            # https://www.alro.com/divsteel/metals_gridpt.aspx?gp=0124&gpn=15-5&Mat=STAINLESS%20STEEL&Type=Bars
            "Steel 15-5": {
                "Yield strength": 1.29e9,  # Pa
                "Ultimate strength": 1.44e9,  # Pa
                "Youngs modulus": 196e9,  # Pa
                "Density": 7800,  # kg/m^3
                "Thermal expansion coefficient": 11e-6  # 1/K
            },
            # https://www.aksteel.nl/files/downloads/172888_armco_17-4_ph_pdb_euro_final_secured_89.pdf
            "Steel 17-4": {
                "Yield strength": 1.267e9,  # Pa
                "Ultimate strength": 1.379e9,  # Pa
                "Youngs modulus": 197e9,  # Pa
                "Density": 7750,  # kg/m^3
                "Thermal expansion coefficient": 11e-6  # 1/K
            }
        }
VEHICLE_MATERIAL = MATERIALS["Aluminium-2024-T4"]
FASTENER_MATERIAL = MATERIALS["Steel 15-5"]

SAFETY_FACTOR = 1.25

# Predefined dimensions
THREAD_PITCH = 0.5e-3  # m; for M3
FORK_CLEARANCE = 2e-3  # m
FASTENER_DIAMETER = 3e-3  # m; M3
FIX_FASTENER_DIAMETER = True

# Initial values
INITIAL_DESIGN_PARAMETERS = {
    "t1": 0,
    "t2": 1e-5,
    "w": 5e-3,
    "x_width": 1e-3,
    "D1": 1e-3,
    "D2": np.zeros(2),  # Pre-defined at this point
    "xs": np.zeros(2),
    "zs": np.zeros(2),
    "lug x offset": FORK_CLEARANCE,
    "min fastener outer": 0,  # Not a design dimension, but data that needs to be carried
}

#INITIAL_FORCES = {
#    "Fy": 642,   # N
#    "Fz": 642,   # N
#    "Fx": 2567,  # N
#    "My": 0,     # Nm
#}

INITIAL_FORCES = {
    "Fz": 642,      # N
    "Fx": 642,      # N
    #"Fx min": 214,  # N
    "Fy": 2567,     # N
    "Mz": 0,        # Nm
}

# External forces
MAX_FORCE_Y = 2139  # N
MAX_FORCE_IN_PLANE = 642  # N
MAX_FORCE_X_BASE = 856  # N
MAX_TOTAL_MOMENT = 0.33 * 2567  # Nm

# Non-lug dimensions
DIFF = 0.18  # m
VEHICLE_WALL_THICKNESS = 21.6e-3  # m

DELTA_T = -47.5  # K

# Precision control
DESIGN_ITERATIONS = 5  # Iterations on the whole design. Chosen "by eye"
THICKNESS_STEPS = 200  # Number of thickness steps
WIDTH_STEPS = 20  # Number of width steps per iteration
WIDTH_ITERATIONS = 5  # Number of width precision iterations


def forces_update(t1):
    """
    :param t1: thickness of the flange
    :return: new z-moment to be carried by the lug, new y force to be carried by the lug.
    Uses deflection to solve for statically indeterminate problem. For this, we need the exact dimensions
    """
    force_x = MAX_TOTAL_MOMENT / (2 * DIFF + t1*t1 / DIFF / 6)
    moment = MAX_TOTAL_MOMENT * t1*t1 / (24 * DIFF * DIFF + 12 * t1*t1)
    return moment, force_x + MAX_FORCE_X_BASE/4


def iteration(material, lug_forces, design_parameters):
    """
    :param material: lug material
    :param lug_forces: forces acting on one lug as a dictionary. WILL BE MANIPULATED !
    :param design_parameters: design parameters as dictionary. WILL BE MANIPULATED !
    :return: mass of the current design in grams (as reference, real effect is the change on dictionaries)
    """
    # Lug dimensions
    allowable_stress = material["Yield strength"] / SAFETY_FACTOR
    best, d1 = flange_design.calculate_flange_dimensions(allowable_stress, lug_forces, design_parameters,
                                                   THICKNESS_STEPS, WIDTH_ITERATIONS, WIDTH_STEPS)

    t1, w, _ = best

    # Define the backplate dimensions
    fastener_inner_diameter = FASTENER_DIAMETER
    fastener_outer_diameter = w
    x_width = max(6 * fastener_inner_diameter, 3 * fastener_inner_diameter + fastener_outer_diameter + t1)
    t2 = design_parameters["t2"]

    # 'Manual' placement of fasteners
    xs = np.array([w/2, x_width - w/2])
    zs = np.array([w/2, w/2])
    d2 = np.ones(2) * fastener_inner_diameter

    effective_flange_area = (w * (w/2 + FORK_CLEARANCE) + w*w/8 * np.pi - d1*d1/4) * t1 / t2

    # Forces on each fastener
    Fx, Fy, Fz = fastener_forces.fastener_forces(w, x_width, lug_forces, xs, zs, d2, effective_flange_area)

    design_parameters["w"] = w
    design_parameters["D2"] = d2
    F_thermal = thermal.thermally_induced_loads(
        DELTA_T, FASTENER_MATERIAL, VEHICLE_MATERIAL, material, design_parameters, VEHICLE_WALL_THICKNESS, THREAD_PITCH
    )

    design_f_thermal = abs(F_thermal)
    Fx += np.ones(2) * design_f_thermal

    if FIX_FASTENER_DIAMETER:
        w = max(w, 3 * fastener_inner_diameter)
        x_width = max(6 * fastener_inner_diameter, 3 * fastener_inner_diameter + fastener_outer_diameter + t1)
    else:
        fastener_inner_diameter = w / 3

    # Pull-through and bearing check check
    pt_margin_lug, pt_margin_vehicle = stress_checks.bearing_check(fastener_inner_diameter, fastener_outer_diameter,
        t2, VEHICLE_WALL_THICKNESS, Fx, Fy, Fz, allowable_stress, VEHICLE_MATERIAL["Yield strength"] * SAFETY_FACTOR)

    min_fastener_outer = 0
    if pt_margin_vehicle < 0:
        print("Vehicle wall fails", pt_margin_vehicle)
        # Minimum fastener head diameter (only parameter to adjust when vehicle wall fails)
        min_fastener_outer = fastener_outer_diameter / (pt_margin_vehicle + 1) - fastener_inner_diameter
        w = fastener_outer_diameter = max(w, min_fastener_outer)
        x_width = max(6 * fastener_inner_diameter, 3 * fastener_inner_diameter + fastener_outer_diameter + t1)

    if pt_margin_lug < 0:
        t2 *= 1 / (pt_margin_lug + 1)

    # Shear stress
    backplate_margin = stress_checks.backplate_check(lug_forces, t1, t2, w, w/2 + FORK_CLEARANCE, allowable_stress)

    if backplate_margin < 0:
        # Not elaborated in the report, since neither asked of, nor influences the design
        print("Backplate becomes an issue")
        t2 *= 1 / (backplate_margin + 1)

    d2 = np.ones(2) * fastener_inner_diameter

    # Update design parameters
    design_parameters["t2"] = t2
    design_parameters["t1"] = t1
    design_parameters["w"] = w
    design_parameters["x_width"] = x_width
    design_parameters["D1"] = d1
    design_parameters["D2"] = d2
    design_parameters["xs"] = xs
    design_parameters["zs"] = zs
    design_parameters["min fastener outer"] = min_fastener_outer

    m = flange_design.mass(design_parameters, material["Density"], False) * 1000  # converted to grams

    # Update Forces
    lug_forces["Mz"], lug_forces["Fy"] = forces_update(t1)

    return m


def run_for_material(material_name):
    """
    For a given material, will run through the iterations and print the mass and design parameters of the final result
    """
    material = MATERIALS[material_name]
    local_forces = INITIAL_FORCES.copy()
    design_parameters = INITIAL_DESIGN_PARAMETERS.copy()
    mass = 0
    for i in range(DESIGN_ITERATIONS):
        mass = iteration(material, local_forces, design_parameters)
    #mass = flange_design.mass(design_parameters, material["Density"], True, True)*1000
    print(material_name, f"mass = {mass:5.2f} g")
    display_params(design_parameters)


def display_params(design_parameters):
    """
    Prints design_parameters in a way easy to the eye
    """
    print(f"""
        t1: {design_parameters['t1'] * 1000:5.2f} mm
        t2: {design_parameters['t2'] * 1000:5.2f} mm
        w : {design_parameters['w'] * 1000:5.2f} mm
        width in x-axis: {design_parameters['x_width'] * 1000:5.2f} mm
        D1: {design_parameters['D1'] * 1000:5.2f} mm
        D2: {design_parameters['D2'][0] * 1000:5.2f} mm
        fasteners at:
            ({design_parameters['xs'][0] * 1000:5.2f}; {design_parameters['zs'][0] * 1000:5.2f}) mm
            ({design_parameters['xs'][1] * 1000:5.2f}; {design_parameters['zs'][1] * 1000:5.2f}) mm
    """[1:])


if __name__ == '__main__':
    #iteration(VEHICLE_MATERIAL, INITIAL_FORCES, INITIAL_DESIGN_PARAMETERS)
    for mat_name in MATERIALS:
        run_for_material(mat_name)
