import numpy as np
"""
 coordinate system is defined as the python coordinate system.
 this means z is positive downwards
 and x positive to the right with (0,0) top left.
"""


def generate_grid(height, width, D, nfastenersx=10, nfastenersz=3):
    """
    Unused
    """
    # define diameter on top row
    Ds = []
    # define number of fasteners in x and z direction.
    # increase in diameters going down in rows, type 0.0 if they are always the same
    # and -x if the diameter decreases when z increases.
    for z in range(nfastenersz):
        for x in range(nfastenersx):
            Ds.append(D)
        D = D + 0.1

    # offset from sides of panel is 1.5 diameter from center of fastener.
    # Functional change: In practice, there is not always enough space to do this the old way -- Nick
    xstart = min(1.5 * Ds[-1], width / 2)
    zstart = min(1.5 * Ds[-1], height / 2)
    # define dimensions panel.
    nfasteners = nfastenersx * nfastenersz
    # equal spacing in both x and z direction.
    xspacing = (width - 2 * xstart) / (nfastenersx - 1) if nfastenersx > 1 else 1e10
    zspacing = (height - 2 * zstart) / (nfastenersz - 1) if nfastenersz > 1 else 1e10
    xend = (nfastenersx - 1) * xspacing + xstart
    zend = (nfastenersz - 1) * zspacing + zstart
    x = np.arange(xstart, xend + xspacing, xspacing)
    z = np.arange(zstart, zend + zspacing, zspacing)

    # Change: to use later -- Nick
    xs = np.zeros(nfasteners)
    zs = np.zeros(nfasteners)

    # the coordinates are defined as the rows consecutively. you can see it if
    # you print the xs and zs list.
    k = 0
    for j in range(nfastenersz):
        for i in range(nfastenersx):
            xs[k] = x[i]
            zs[k] = z[j]

    return xs, zs, np.array(Ds)


def fastener_forces(width, length, forces, xs, zs, Ds, a_effective_lug):
    """
    :param width: w
    :param length: x-width of the lug backplate
    :param forces: current lug-forces dictionary
    :param xs: x-positions for each fastener hole, as array
    :param zs: z-positions for each fastener hole, as array
    :param Ds: diameters, as array
    :param a_effective_lug:
    :return: x-, y-, z-forces on each fastener, as array
    Calculated with asymetric design in mind. Could be simplified to 8 lines
    """

    Fx, Fy, Fz, moment = forces["Fx"], forces["Fy"], forces["Fz"], forces["Mz"]

    nfasteners = len(xs)

    area = np.pi * Ds * Ds / 4

    # this is the standard cg location formula, divided into two steps.
    numeratorx = sum(area * xs)
    numeratorz = sum(area * zs)

    # Functional edit: the mass is the one of the whole lug. Also have to consider the effective area of the lug -- Nick
    A_tot = width * length + a_effective_lug
    x_cg = (A_tot * length/2 - numeratorx) / (A_tot - sum(area))
    z_cg = (A_tot * width/2 - numeratorz) / (A_tot - sum(area))

    centrelugx = length / 2
    centrelugz = width / 2
    # determine the moment

    # 2 determine the moment around the centre of the lug (acw+)
    dx = centrelugx - x_cg
    dz = centrelugz - z_cg
    My = Fx * dz - Fz * dx

    # Part 3 is determining the force induced by this moment on each fastener
    # 1 determine the distance of the individual fasteners to the cog

    diff_x, diff_z = xs - np.ones(nfasteners) * x_cg, zs - np.ones(nfasteners) * z_cg
    r = np.sqrt(diff_x*diff_x + diff_z*diff_z)

    # 2 determine the sum of the area's multiplied with the distance squared
    Ar2 = sum(area * r*r)

    # 3 determine the net force on each fastener by the moment
    F = My * area * r / Ar2

    angle = np.arctan(diff_z * np.reciprocal(diff_x))
    Fxi = np.cos(angle) * F + np.ones(nfasteners) * Fx / nfasteners
    Fzi = np.sin(angle) * F + np.ones(nfasteners) * Fz / nfasteners

    div = sum(moment * xs*xs)
    if div != 0:
        Fyi = moment * xs / div + np.ones(nfasteners) * Fy / nfasteners
    else:
        Fyi = np.ones(nfasteners) * Fy / nfasteners

    return np.array(Fxi), Fyi, np.array(Fzi)
