import numpy as np
# import matplotlib.pyplot as plt

from math import pi, sqrt
import json

H2_PRESSURE_MIN = .07703e5
H2_PRESSURE_MAX = 13.301e5


def calc_thickness(material, launch, masses, dimensions, internal_pressure):
    allowable_stress = material["Yield strength"] / 1.25
    density = material["Density"]
    radius = dimensions["Radius"]
    length = masses["H2 Fuel"] / (dimensions["Fuel density"] * pi * radius**2) + 2/3 * radius

    acc = sqrt(launch["Max Acceleration total"]**2 + launch["Lateral Acceleration"]**2)
    hydrostatic_pressure = acc * dimensions["Fuel density"] * length
    p = np.linspace(hydrostatic_pressure + H2_PRESSURE_MIN, hydrostatic_pressure + H2_PRESSURE_MAX, 100)
    thicknesses = p * radius / allowable_stress
    m = mass(material, radius, length, thicknesses)
    req_pressure = acc * thicknesses * density * 2 * (length / radius - 1) + np.ones(100) * acc * masses["Communication"] / (pi * radius**2)

    # This are the actual values
    p = hydrostatic_pressure + internal_pressure
    t1 = p * ((length - 2 * radius) * 2 * radius + pi * radius ** 2) / (allowable_stress * (2 * (length - 2 * radius) + 2 * pi * radius))
    m1 = mass(material, radius, length, t1)
    comp_p = 9.81 * (m1 + masses["Communication"]) / (pi * radius**2)
    #print(p1)

    #plt.plot(p * 1e-5, weights)
    #plt.plot(p * 1e-5, req_pressure * np.reciprocal(p))
    #plt.show()
    return radius, length, t1


def mass(material, radius, length, thickness):
    return (4 * pi * radius * radius + 2 * pi * radius * length) * thickness * material["Density"]


if __name__ == '__main__':
    with open("../current_parameters.json", "r") as f:
        data = json.load(f)

    launch = data["Launch"]
    sc_masses = data["Spacecraft masses"]
    dimensions = data["FuelTank"]

    print(calc_thickness(data["Materials"]["Cryo Ti-6AI-4V"], launch, sc_masses, dimensions, dimensions["Internal pressure"]))
