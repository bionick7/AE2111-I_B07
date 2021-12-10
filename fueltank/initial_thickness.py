import numpy as np
# import matplotlib.pyplot as plt

from math import pi, sqrt
import json

"""
RADIUS = 2.2
FUEL_MASS = 5313.170
M_COMMS = 121.9

ACCELERATION = 6.32455532 * 9.81
H2_DENSITY = 70.850  # @ 70.017 K, 1 bar
TI_STRENGTH = 1590e6 / 1.25
TI_DENSITY = 4430

H2_PRESSURE = 101325
"""

H2_PRESSURE_MIN = .07703e5
H2_PRESSURE_MAX = 13.301e5


def calc_thickness(material, launch, masses, dimensions):
    allowable_stress = material["Yield strength"] / 1.25
    density = material["Density"]
    radius = dimensions["Radius"]
    length = masses["H2 Fuel"] / (dimensions["Fuel density"] * pi * radius**2) + 2/3 * radius

    acc = sqrt(launch["Max Acceleration total"]**2 + launch["Lateral Acceleration"]**2)
    hydrostatic_pressure = acc * dimensions["Fuel density"] * length
    p = np.linspace(hydrostatic_pressure + H2_PRESSURE_MIN, hydrostatic_pressure + H2_PRESSURE_MAX, 100)
    thicknesses = p * radius / allowable_stress
    weight = (4 * pi * radius * radius + 2 * pi * radius * length) * thicknesses * density
    req_pressure = acc * thicknesses * density * 2 * (length / radius - 1) + np.ones(100) * acc * masses["Communication"] / (pi * radius**2)

    # This are the actual values
    p = hydrostatic_pressure + dimensions["Internal pressure"]
    t1 = p * radius / allowable_stress
    w1 = (4 * pi * radius * radius + 2 * pi * radius * length) * t1 * density
    comp_p = 9.81 * (w1 + masses["Communication"]) / (pi * radius**2)
    #print(p1)

    #plt.plot(p * 1e-5, weights)
    #plt.plot(p * 1e-5, req_pressure * np.reciprocal(p))
    #plt.show()
    return w1, t1


if __name__ == '__main__':
    with open("../current_parameters.json", "r") as f:
        data = json.load(f)

    launch = data["Launch"]
    sc_masses = data["Spacecraft masses"]
    dimensions = data["FuelTank"]

    print(calc_thickness(data["Materials"]["Cryo Ti-6AI-4V"], launch, sc_masses, dimensions))
