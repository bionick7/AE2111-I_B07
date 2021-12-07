import numpy as np
import matplotlib.pyplot as plt

from math import pi

RADIUS = 2.2
FUEL_MASS = 5313.170
M_COMMS = 121.9

ACCELERATION = 6.32455532 * 9.81
H2_DENSITY = 70.850  # @ 70.017 K, 1 bar
TI_STRENGTH = 1590e6 / 1.25
TI_DENSITY = 4430

H2_PRESSURE_MIN = .07703e5
H2_PRESSURE = 101325
H2_PRESSURE_MAX = 13.301e5


def get_weight_thickness(p, length):
    thickness = p * RADIUS / TI_STRENGTH
    weight = (4 * pi * RADIUS * RADIUS + 2 * pi * RADIUS * length) * thickness * TI_DENSITY
    return weight, thickness


length = FUEL_MASS / (H2_DENSITY * pi * RADIUS**2) + 2/3 * RADIUS
print(length)

hydrostatic_pressure = ACCELERATION * H2_DENSITY * length
p = np.linspace(hydrostatic_pressure + H2_PRESSURE_MIN, hydrostatic_pressure + H2_PRESSURE_MAX, 100)
weights, thicknesses = get_weight_thickness(p, length)
req_pressure = ACCELERATION * thicknesses * TI_DENSITY * 2 * (length / RADIUS - 1) + np.ones(100) * ACCELERATION * M_COMMS / (pi * RADIUS**2)

w1, t1 = get_weight_thickness(hydrostatic_pressure + H2_PRESSURE, length)
p1 = 9.81 * (w1 + M_COMMS) / (pi * RADIUS**2)
print(w1, t1*1000)
print(p1)

#plt.plot(p * 1e-5, weights)
#plt.plot(p * 1e-5, req_pressure * np.reciprocal(p))
#plt.show()

