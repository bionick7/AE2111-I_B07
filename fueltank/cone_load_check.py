import math
import json
import numpy as np
import matplotlib.pyplot as plt

R_top_outer = 1.80641
cg_arm = 3.444  # todo does this change?
#m = 5582.9
g = 9.81
#yield_strength = 500 * 10**6

DIVISIONS = 100

ANGLE = math.radians(34.805)


def calculate_thicknesses(material, Ftrans, moment_top, F_max_compress):
    Torque = 0.04

    yield_strength = material["Yield strength"]

    t_list = np.zeros(DIVISIONS)
    x_list = np.zeros(DIVISIONS)
    for i, x in enumerate(np.linspace(0, 0.94426, DIVISIONS)):
        r = R_top_outer - math.tan(ANGLE) * x
        t = 0.01e-3

        stress_total = math.inf
        while True:
            if stress_total <= yield_strength:
                break
            t *= 1.005
            sigma = (F_max_compress/(2 * math.pi * r) + moment_top/(math.pi * r**2))/t
            true_compression, shear_compression = sigma * math.cos(ANGLE), sigma * math.sin(math.radians(34.805))
            shear = (2 * Ftrans)/(math.pi * r * t) + Torque/(2 * math.pi * t * r**2)

            stress_total = math.sqrt(true_compression**2 + 3 * shear_compression**2 + 3 * shear**2)

        t_list[i] = t
        x_list[i] = x
        #print(stress_total/1e6, '[MPa]')
        #print(c)

    return x_list, t_list


if __name__ == '__main__':

    with open("../current_parameters.json", "r") as f:
        data = json.load(f)

    x_list, t_list = calculate_thicknesses(data["Materials"]["Ti-6AI-4V"], 5582.9)

    plt.plot(x_list, t_list)
    plt.ylabel('Thickness [mm]')
    plt.xlabel('Distance from top [m]')


    plt.show()

    #I_x = 0.25 * math.pi * y**4
    #A_top = math.pi * (R_top_outer**2 - (R_top_outer - t)**2)
