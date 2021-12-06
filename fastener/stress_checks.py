from math import pi
import numpy as np


def bearing_check(inner_diameter, outer_diameter, t2, vehicle_wall_thickness, Fx, Fy, Fz, allowable_lug, allowable_vehicle):
    """
    :param inner_diameter: inner diameter of the bolt
    :param outer_diameter: outer diameter of the bolt
    :param t2: wall thickness of the lug
    :param vehicle_wall_thickness: thickness of the wall, the lug is attached to
    :param Fx: Force in the x-axis per fastener
    :param Fy: Similar for y
    :param Fz: Similar for z
    :return: The margins for the lug and the vehicle (in that order)
    """

    shear_stress_lug = Fx * (1 / (t2 * inner_diameter * pi))
    shear_stress_vehicle = Fx * (1 / (vehicle_wall_thickness * inner_diameter * pi))
    normal_stress = Fx / (np.pi * (outer_diameter - inner_diameter)**2)
    shear_stress = np.sqrt(Fy * Fy + Fz * Fz)
    bearing_stresses_lug = shear_stress * np.reciprocal(inner_diameter) / t2
    bearing_stresses_vehicle = shear_stress * np.reciprocal(inner_diameter) / vehicle_wall_thickness

    tension_yield_stress_lug = max(np.sqrt(0.5 * (normal_stress**2 + bearing_stresses_lug**2) + 3 * shear_stress_lug**2))
    tension_yield_stress_vehicle = max(np.sqrt(0.5 * (normal_stress**2 + bearing_stresses_vehicle**2) + 3 * shear_stress_vehicle**2))
    return allowable_lug / tension_yield_stress_lug - 1, allowable_vehicle / tension_yield_stress_vehicle - 1


def backplate_check(lug_forces, t1, t2, w, moment_arm, allowable):
    area = t2 * w * 2

    f_z = lug_forces["Fz"]
    f_x = lug_forces["Fx"]
    f_y = lug_forces["Fy"] + lug_forces["Fx"] * moment_arm / t1

    shear_stress_lug_1 = f_z / area      # Both can be either y or z, depending on the lug
    shear_stress_lug_2 = f_x / area  # Both can be either y or z, depending on the lug
    normal_stress_lug = f_y / area       # Always x
    stress = np.sqrt(0.5 * normal_stress_lug**2 + 3 * shear_stress_lug_1**2 + 3 * shear_stress_lug_2**2)
    return allowable / stress - 1


def fastener_check(Fx, Fy, Fz, inner_diameter, allowable):
    """
    :param Fx: forces in the x-direction for each fastener, as an array
    :param Fy: forces in the y-direction for each fastener, as an array
    :param Fz: forces in the z-direction for each fastener, as an array
    :param inner_diameter: fastener / hole diameter
    :param allowable: maximal allowable stress for fastener
    :return: The stress margin for the dastener
    Perform a simple check of the material cylinder, does not take into account thread damage etc.
    Unused
    """
    fastener_area = inner_diameter * inner_diameter * np.pi / 4
    fastener_normal = Fx / fastener_area
    fastener_shear = np.sqrt(Fy*Fy + Fz*Fz) / fastener_area
    fastener_stress = max(np.sqrt(0.5 * fastener_normal**2 + 3 * fastener_shear**2))
    return allowable / fastener_stress - 1

