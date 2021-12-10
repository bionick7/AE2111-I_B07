# For each material:
#   Calculate initial dimensions
#   Verify initial buckling
#   Calculate cone dimensions
#   Verify cone buckling

import json
import matplotlib.pyplot as plt
import numpy as np

import initial_thickness
import cone_load_check
import buckling


ANGLE = 34.805


def cone_volume(xs, ts):
    L = 1.15
    dx = L / len(xs)
    terms = 2 * np.pi * (L + np.sin(ANGLE) * xs) * ts
    return np.sum(terms * dx) + (terms[0] + terms[-1]) / 2 * dx  # Trapezoidal rule https://en.wikipedia.org/wiki/Trapezoidal_rule


def main():
    with open("../current_parameters.json", "r") as f:
        data = json.load(f)

    launch = data["Launch"]
    sc_masses = data["Spacecraft masses"]
    dimensions = data["FuelTank"]
    ft_material = data["Materials"]["Cryo Ti-6AI-4V"]
    connector_material = data["Materials"]["Ti-6AI-4V"]

    internal_pressure = 101325
    prev_t = 0
    R, L, t = 0, 0, 1
    mass = 0
    while abs(t - prev_t) > 1e-5:
        prev_t = t
        R, L, t = initial_thickness.calc_thickness(ft_material, launch, sc_masses, dimensions, internal_pressure)
        mass = initial_thickness.mass(ft_material, R, L, t)

        total_mass = mass + sc_masses["Communication"] + sc_masses["H2 Fuel"]
        load = launch["Max Acceleration total"] * total_mass * 9.81
        moment = launch["Lateral Acceleration"] * total_mass * 9.81 * 3.444

        load_factor = 1
        delta_p = 0
        while load_factor > 0.6667:
            buckling_force, buckling_moment = buckling.cone_buckling_bending(ft_material, ANGLE, t, R, delta_p)
            load_factor = load / buckling_force + moment / buckling_moment
            delta_p += 10

        internal_pressure = 101325 + delta_p

    total_mass = mass + sc_masses["Communication"] + sc_masses["H2 Fuel"]
    F_max_adcs = 88 * (2 ** 0.5)
    compressive = launch["Max Acceleration total"] * total_mass * 9.81 + F_max_adcs  # todo It won't fire ADCS during launch ?!
    transverse = launch["Lateral Acceleration"] * total_mass * 9.81 # ((8) ** 0.5) * g * m
    moment = transverse * 3.444

    xs, ts = cone_load_check.calculate_thicknesses(connector_material, transverse, moment, compressive)

    load_factor = 10
    t_min = min(ts)
    while load_factor > .6667:
        buckling_force, buckling_moment = buckling.cone_buckling_bending(connector_material, ANGLE, t_min, R, 0)
        load_factor = compressive / buckling_force + moment / buckling_moment
        t_min += 1e-4
    print(t_min)
    ts_2 = np.maximum(t_min * np.ones(len(ts)), ts)

    print(cone_volume(xs, ts) * connector_material["Density"])
    print(cone_volume(xs, ts_2) * connector_material["Density"])

    plt.plot(xs, ts * 1000)
    plt.plot(xs, ts_2 * 1000)
    plt.show()


if __name__ == '__main__':
    main()
