import math
import numpy as np
import matplotlib.pyplot as plt

R_top_outer = 1.80641
cg_arm = 3.444
m = 5582.9
g = 9.81
yield_strength = 500 * 10**6

Ftrans = ((8)**0.5) * g * m
moment_top = cg_arm * Ftrans

F_max_adcs = 88 * (2**0.5)
F_max_compress = ((m * g * 6) + F_max_adcs)
Torque = 0.04

print(Ftrans)
t_list = []
x_list = []


for x in np.linspace(0, 0.94426, 20):
    r = R_top_outer - math.tan(math.radians(34.805))* x 
    c = 0.01

    stress_total = math.inf
    while True:
        if stress_total <= yield_strength:
            break
        c = c*1.005
        t = c/1000
        sigma = (F_max_compress/(2 * math.pi * r) + moment_top/(math.pi * r**2))/t
        stress_total = ((sigma * math.cos(math.radians(34.805)))**2 + 3 * (math.sin(math.radians(34.805)) * sigma)**2 + 3 * ((2 * Ftrans)/(math.pi * r * t) + Torque/(2 * math.pi * t * r**2))**2)**0.5

    t_list.append(c)
    x_list.append(x)
    print(stress_total/1e6, '[MPa]')
    print(c)

plt.plot(x_list, t_list)
plt.ylabel('Thickness [mm]')
plt.xlabel('Distance from top [m]')


plt.show()

#I_x = 0.25 * math.pi * y**4
#A_top = math.pi * (R_top_outer**2 - (R_top_outer - t)**2)
