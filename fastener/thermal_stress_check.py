#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Nov 24 19:21:41 2021

@author: laurahoitingh
"""

from math import pi


def force_ratio(E_backplate, E_spacecraft, E_fastener, design_parameters, spacecraft_wall_thickness, thread_pitch):
    """
    :params E_...: Youngs modulus of the lug material
    :return: The force ratio for the spacecraft and the spacecraft wall
    Used by thermally_induced_loads
    """
    t2 = design_parameters["t2"]
    Dfo = design_parameters["w"]
    Dfi = design_parameters["D2"][0]

    Lhsub = 0.4 * Dfi  # assuming a cylindrical head
    Lengsub = 0.4 * Dfi  # assuming a nut
    Lnsub = 0.4 * Dfi
    Lshank = 0.02  # Part of the bolt which has no threading and thus uses full diameter
    Lthreaded = t2 + spacecraft_wall_thickness - Lshank
    L = [Lhsub, Lengsub, Lnsub, Lshank, Lthreaded]

    Ah = pi * Dfo**2 / 4
    Anom = pi * Dfi**2 / 4
    Athreaded = (pi / 4) * (Dfi - 0.938194 * thread_pitch)**2
    A = [Ah, Anom, Anom, Anom, Athreaded]

    da2 = (4 * t2) / ((E_backplate * pi) * (Dfo**2 - Dfi**2))
    da3 = (4 * spacecraft_wall_thickness) / ((E_spacecraft * pi) * (Dfo**2 - Dfi**2))
    LowerA = [x / y for x, y in zip(L, A)]
    sumLowerA = sum(LowerA)
    db = (1/E_fastener) * sumLowerA
    phi_lug = da2 / (da2 + db)
    phi_sc_wall = da3 / (da3 + db)
    return phi_lug, phi_sc_wall


def thermally_induced_loads(delta_T, fastener_material, spacecraft_material, lug_material, design_parameters,
                            spacecraft_wall_thickness, thread_pitch):
    """
    :return: The force ratio for the spacecraft and the spacecraft wall
    """
    phi_lug, phi_sc_wall = force_ratio(lug_material["Youngs modulus"], spacecraft_material["Youngs modulus"],
                        spacecraft_material["Youngs modulus"], design_parameters, spacecraft_wall_thickness, thread_pitch)
    # coefficients of thermal expansion and Young's modulus for fastener
    alpha_b, E_b = fastener_material["Thermal expansion coefficient"], fastener_material["Youngs modulus"]
    # coefficient of thermal expansion of spacecraft wall material
    alpha_w = spacecraft_material["Thermal expansion coefficient"]

    # coefficients of thermal expansion for clamped parts
    alpha_c = (lug_material["Thermal expansion coefficient"] * design_parameters["t2"] + alpha_w * spacecraft_wall_thickness) \
              / (design_parameters["t2"] + spacecraft_wall_thickness)
    
    # stiffness area of fastener:
    A_sm = 0.25 * pi * (design_parameters["D2"][0] - 1.22687 * thread_pitch)**2
    
    # thermal induced load:
    F_deltaT = (alpha_c - alpha_b) * delta_T * E_b * A_sm * (1 - phi_sc_wall)  # Use phi lug or phi S/C ???
    
    return F_deltaT
