import math
import json

# All values labeled #PLACEHOLDER are still incorrect values. They are only there to allow the program to be test-run, but ARE NOT CORRECT for the S/C.

G = 9.81


def Euler_column(material, dimensions):
    """
    :param material: Dict containing material properties
    :param dimensions: Dict containing the dimensions of the tank
    
    Returns the critical euler column buckling stress, derived from the given dimensions and material properties.
    WARNING: Does NOT yet include any safety factors.
    """
    # Extract the parameters from the dicts
    l_tank = dimensions["Length"]
    radius = dimensions["Radius"]
    thickness = dimensions["Wall thickness"]
    l_cylinder = l_tank - 2 * radius  # Calculate the length of the cylindrical section of the cylinder.

    E = material["Modulus"]
    # Calculate the material properties.
    I = 1/4 * math.pi * (radius ** 4)
    A_encl = math.pi * radius * radius  # "Enclosed" area, aka the area of the vertical projection of the tank.
    sigma_crit = (math.pi * math.pi * E * I) / (A_encl * l_cylinder * l_cylinder)
    return sigma_crit


def Shell_buckling(material, dimensions):
    """
    :param material: Dict containing material properties
    :param dimensions: Dict containing the dimensions of the tank
    :param requirements: Dict containing the requirements for the tank wall
    
    Returns the critical shell buckling stress, derived from the given dimensions and material properties".
    WARNING: Does NOT yet include any safety factors.
    """
    l_tank = dimensions["Length"]
    radius = dimensions["Radius"]
    thickness = dimensions["Wall thickness"]

    delta_p = dimensions["Internal pressure"] - 101325

    E = material["Modulus"]
    poisson = material["Poisson"]
    l_cylinder = l_tank - 2 * radius #Calculate the length of the cylindrical section of the cylinder.

    lambd = l_cylinder**2 * math.sqrt(12 * (1 - poisson**2)) / (math.pi**2 * radius * thickness)
    k = 2 * lambd #Love this. Minimum value of y = lambd + a / lambd, where a == lambd^2

    Q = delta_p / E * (radius / thickness) ** 2

    sigma_crit = (1.983 - 0.983 * math.exp(-23.14 * Q)) * k * math.pi ** 2 * E / (12 * (1 - poisson**2)) * (thickness / l_tank)**2
    return sigma_crit


def Max_cylinder_wall_stress(launch, sc_masses, dimensions):
    """
    Returns the maximum stress (both compressive and tensile) in the cylindrical part of the wall.
    """
    long_g_max = launch["Max Acceleration early launch"] * G
    long_g_min = launch["Min Acceleration"] * G
    lateral_g = launch["Lateral Acceleration"]
    
    mass_ant = sc_masses["Antenna"]
    mass_tank = dimensions["Mass"]
    mass_prop = sc_masses["H2 Fuel"]

    radius = dimensions["Radius"]
    thickness = dimensions["Wall thickness"]
    l_tank = dimensions["Length"]
    d_ant = dimensions["Antenna Distance"] + l_tank - radius  # The distance from the bottom of the cylindrical structure.
    l_cylinder = l_tank - 2 * radius

    A_wall = 2 * math.pi * radius * thickness  # Thin wall approximation of total wall area
    A_encl = math.pi * radius * radius
    delta_p = dimensions["Internal pressure"] - 101325  # That's just it
    I_xx = math.pi * radius**3 * thickness  # Area moment of inertia
    
    # Stresses:
    # Longitudinal acceleration: sigma = m * a / A
    # Lateral acceleration: sigma = m * l * a * r / I_xx
    # Pressure: sigma = delta_p * A_encl / A
    # Propellant mass acceleration (only for the deceleration / tension case, as then the force acts on the top dome.
    
    # Compression:
    sigma_comp_long = - long_g_max * (mass_ant + mass_tank) / A_wall
    sigma_comp_lat = - ((mass_ant * d_ant) + ((mass_tank + mass_prop) * l_cylinder / 2)) * lateral_g * radius / I_xx
    sigma_comp_press = A_encl * delta_p / A_wall
    sigma_comp = sigma_comp_long + sigma_comp_lat + sigma_comp_press

    # Tension
    sigma_tens_long = - long_g_min * (mass_ant + mass_tank + mass_prop) / A_wall
    sigma_tens_lat = - sigma_comp_lat
    sigma_tens_press = (A_encl * dimensions["Internal pressure"] / A_wall)
    sigma_tens = sigma_tens_long + sigma_tens_lat + sigma_tens_press
    
    # Probably just return both.
    return sigma_comp, sigma_tens


if __name__ == '__main__':

    with open("../current_parameters.json", "r") as f:
        data = json.load(f)

    launch = data["Launch"]
    sc_masses = data["Spacecraft masses"]
    dimensions = data["FuelTank"]

    print(Shell_buckling(data["Materials"]["Ti-6AI-4V"], dimensions) / -Max_cylinder_wall_stress(launch, sc_masses, dimensions)[0])
