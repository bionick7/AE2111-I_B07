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

BUCKLING_SF = 1.5

ANGLE = 34.805 * np.pi / 180
ANGLE2 = 41.288 * np.pi / 180


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
    connector_material = data["Materials"]["Aluminium-2024-T4"]
    connector_material_it1 = data["Materials"]["Ti-6AI-4V"]

    ft_mass = design_fueltank(ft_material, launch, sc_masses, dimensions)
    connector_mass = design_connector(connector_material, launch, sc_masses, ft_mass, connector_material_it1)
    total_mass = ft_mass + connector_mass
    print(ft_mass)


def design_fueltank(ft_material, launch, sc_masses, dimensions):
    internal_pressure = 101325
    prev_t = 0
    R, L, t = 0, 0, 1
    mass = 0
    allowable_stress = ft_material["Yield strength"] / 1.25
    while abs(t - prev_t) > 1e-5:
        prev_t = t
        R, L, t = initial_thickness.calc_thickness(ft_material, launch, sc_masses, dimensions, internal_pressure)
        mass = initial_thickness.mass(ft_material, R, L, t)

        total_mass = mass + sc_masses["Communication"] + sc_masses["H2 Fuel"]
        load = launch["Max Acceleration total"] * total_mass * 9.81
        moment = launch["Lateral Acceleration"] * total_mass * 9.81 * 3.444

        load_factor = 1
        delta_p = 0
        while load_factor > 1 / BUCKLING_SF:
            t = (101325 + delta_p) * ((L - 2*R) * 2*R + np.pi * R ** 2) / (allowable_stress * (2 * (L - 2*R) + 2 * np.pi * R))
            buckling_force, buckling_moment = buckling.cone_buckling_bending(ft_material, ANGLE, t, R, delta_p)
            load_factor = load / buckling_force + moment / buckling_moment
            delta_p += 1

        internal_pressure = 101325 + delta_p

    print(t)

    return mass


def design_connector(connector_material, launch, sc_masses, ft_mass, it1_material):
    total_mass = ft_mass + sc_masses["Communication"] + sc_masses["H2 Fuel"]
    F_max_adcs = 88 * (2 ** 0.5)
    compressive = launch["Max Acceleration total"] * total_mass * 9.81 + F_max_adcs  # todo It won't fire ADCS during launch ?!
    transverse = launch["Lateral Acceleration"] * total_mass * 9.81 # ((8) ** 0.5) * g * m
    moment = transverse * 3.444

    xs, ts = cone_load_check.calculate_thicknesses(connector_material, transverse, moment, compressive)
    _, ts_it1 = cone_load_check.calculate_thicknesses(it1_material, transverse, moment, compressive)

    radius = np.ones(len(xs)) * 1.80641 - np.sin(ANGLE) * xs
    buckling_force, buckling_moment = buckling.cone_buckling_bending(connector_material, ANGLE, ts, radius, 0)
    load_factor = compressive / buckling_force + moment / buckling_moment

    ts2 = ts * np.sqrt(load_factor * 1.5)
    #print(buckling.cone_buckling_bending(connector_material, ANGLE, ts2, radius, 0)[0])

    #print(ts2)

    #print("t =", t_min * 1000, "mm")
    ts_tot = np.maximum(ts, ts2)


    ts_it1 = np.flip(ts_it1, 0)
    ts2 = np.flip(ts2, 0)

    print("m1 =", cone_volume(xs, ts_it1) * it1_material["Density"], "kg")
    print("m2 =", cone_volume(xs, ts2) * connector_material["Density"], "kg\n")

    plt.plot(xs / np.cos(ANGLE), ts_it1 * 1000, label="1$^{st}$ iteration (Titanium)")
    plt.plot(xs / np.cos(ANGLE), ts2 * 1000, label="2$^{nd}$ iteration (Aluminium)")
    plt.xlabel("Location along cone -- bottom to top [m]")
    plt.ylabel("thickness [mm]")
    plt.legend(loc=5)
    plt.axhline(0, c="black", ls="--")
    plt.axvline(0, c="black", ls="-", linewidth=1)
    plt.axvline(1.15, c="black", ls="-", linewidth=1)
    plt.show()

    return cone_volume(xs, ts2) * connector_material["Density"]


if __name__ == '__main__':
    main()
    #S = .16**2 * .17653 * (3.4*2.81) ** .33333 * .01 ** .6667 * (40.18**2) / (0.017 **2) * 482.6/503
    #print(S)