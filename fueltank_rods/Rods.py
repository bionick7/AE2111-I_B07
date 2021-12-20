import math
import numpy as np
from Rod_buckling import *
from Materials import materials
import matplotlib.pyplot as plt


dimensions = {
    #All major dimensions have been found using Catia.
    "Length": 1.150, #Meters. Length of the connecting rod.
    "Minimum Thickness": 1e-4, #Minimum rod thickness in meters
    "CoM Distance": 3.444,
    "Mass": 5582.9,
    "Rod Connection Height": 0.748345, #Height in meters. (Vertical distance between the connection point on the tank, and the connection point on the cylinder).
    "Rod Connection Radius": 1.6531, #Meters
    "Rod Vertical Angle": 40.597, #Degrees
    "Rod Horizontal Angle": 41.185, #Degrees
    "Minimum Diameter": 0.01, #Minimum rod radius in meters
    "Maximum Diameter": 0.5, #Maximum rod radius in meters
    "Diameter Step": 0.0001,
    }

requirements = {
    "Lateral Acceleration": 2, #G
    "Longitudinal Acceleration": -6, #G #Negative as this will require a force positive upwards, which is negative y in the coordinate system used here.
    "Longitudinal Deceleration": 2.5, #G
    "Safety Factor": 1.25, #From Ariane 5 User Manual, required design safety factor for the S/C structure.
    }

def Find_optimal(materials, dimensions, requirements):
    d = dimensions #Let's keep it somewhat readable and compact, shall we?
    optimums = [] #Initialize the list for all optimum values
    forces = Max_force(requirements, dimensions)
    print(forces)
    for material_name in materials: #For each material:
        material = materials[material_name]
        data = [] #Initialize the list of data points for this material
        for diameter in np.arange(d["Minimum Diameter"], d["Maximum Diameter"], d["Diameter Step"]):
            dimensions["Radius"] = diameter / 2
            thickness = find_thickness(material, dimensions, requirements, forces)
            if thickness: #If a valid thickness was found:
                data.append([diameter + thickness, Mass(material, dimensions, thickness), thickness * 1e3]) #Append the data point to the data list
        #Find the most lightweight / optimal design
        masses = [point[1] for point in data]
        thicks = [point[2] for point in data]
        try:
            optimal = sort_by(data, 1)[0] #Sort all items by mass
        except IndexError:
            print(f"No valid configuration found for {material_name}")
        t_opt = optimal[2]
        m_opt = optimal[1]
        d_opt = optimal[0]
        #Print the optimum values
        print(f"Material: {material_name}")
        print(f"Minimum mass achieved at a diameter of {d_opt:.3f} m, with thickness {(1e3 * t_opt):.6f} mm. Mass is {m_opt:.5f} kg.")
        #Add the optimum values to the optimums list
        optimums.append([material_name] + optimal)
        #Add the diameter-mass plot to the plot
        plt.plot([point[0] for point in data], masses, label = material_name)
    #Print the most optimal item
    print(f"Most optimal material: {sort_by(optimums, 2)[0][0]}")
    optimum = sort_by(optimums, 2)[0]
    print(f"Total mass: {6 * optimum[2]:.3f} kg")
    #Finalise the plot
    plt.legend()
    #plt.title("Mass vs rod diameter")
    plt.xlabel("Rod diameter [m]")
    plt.ylabel("Rod mass [kg]")
    #plt.ylabel("Rod wall thickness [mm]")
    plt.show()
    #Return the best item
    return(sort_by(optimums, 2)[0])
            

def Mass(material, dimensions, thickness):
    """
    Returns the mass of the cylinder for the given configuration.
    """
    length = dimensions["Length"]
    radius = dimensions["Radius"]
    return(Area(radius, thickness) * length * material["Density"])


def find_thickness(material, dimensions, requirements, forces):
    sigma_y = material["Yield Strength"]
    radius = dimensions["Radius"]
    l_rod = dimensions["Length"]
    t_min = dimensions["Minimum Thickness"]
    t_max = radius / 10 #The thickness may not exceed 1 / 10 * radius, as that would invalidate the thin-walled assumptions.
    t_max = radius * 2
    #Start the iteration phase:
    t = t_min
    f_comp = abs(forces[0]) #The maximum compressive force, expressed as a positive value
    f_tens = abs(forces[1]) #The maximum tensile force in the rods, also positive
    while True:
        #Select the most critical failure stresses stress (the lowest value is the most critical)
        euler = Euler_column(material, dimensions, t)
        shell = Shell_buckling(material, dimensions, t)
        if radius / t >= 10: #If the thin walled assumption is valid, include Shell buckling. Else exclude it.
            sigma_crit_comp = min([Euler_column(material, dimensions, t), Shell_buckling(material, dimensions, t), sigma_y])
        else:
            sigma_crit_comp = min([Euler_column(material, dimensions, t), sigma_y])
        sigma_crit_tens = min([sigma_y])
        area = Area(radius, t)
        sigma_comp = f_comp / area
        sigma_tens = f_tens / area
        if sigma_comp < sigma_crit_comp and sigma_tens < sigma_crit_tens: #If the actual stress is below the critical stress:
            if False:
                if euler < min(sigma_y, shell):
                    print("Eul")
                elif shell < min(sigma_y, euler):
                    print("She")
                else:
                    print("Yield")
                    #print(euler / sigma_y)
                #print(sigma_comp)
            return(t)
        t *= 1.001 #Increase t by 0.1% if the thickness is too low
        if t > t_max: #If the thin-walled assumption is no longer valid, give up. We arise by the thin walled assumption, and if we have to, so will we perish.
            return(None)

def sin(x):
    """
    Returns the sine of x, in degrees
    """
    return(math.sin(math.radians(x)))

def cos(x):
    """
    Returns the cosine of x, in degrees
    """
    return(math.cos(math.radians(x)))


def Forces(requirements, dimensions):
    """
    Calculates the forces exerted by each of the rods in the "Stewart Platform",
    required to cause a certain F_x, F_z, M_x, ... at the top side of the
    platform. (The central point between the three connection points.) To get
    the forces inside each rod, all forces have to be inverted.
    """
    #Set up the constants:
    theta = dimensions["Rod Vertical Angle"] #Angle between the horizontal plane, and the axis of the rods
    alpha = dimensions["Rod Horizontal Angle"] #Angle between the radial-in and the projection of the rod on the horizontal plane
    r = dimensions["Rod Connection Radius"] #The radius from the center of the tank, to the connection point on the tank
    #Set up the forces
    F_x = requirements["F_x"]
    F_y = requirements["F_y"]
    F_z = requirements["F_z"]
    #Moments around the base of the plate.
    moment_arm = dimensions["CoM Distance"] - dimensions["Rod Connection Height"]
    #Calculate the moments the rods have to exert to counteract the moment caused by any sideways force.
    M_x = -F_z * moment_arm
    M_y = 0
    M_z = F_x * moment_arm

    #Set up the matrices
    #Note: Only the horizontal components.
    F_x_m = cos(theta) * np.array([sin(alpha), -sin(alpha), cos(30 + alpha), cos(30 - alpha), -cos(30 - alpha), -cos(30 + alpha)])
    #Basically F_x_m, but shifted a bit.
    F_z_m = cos(theta) * np.array([cos(alpha), cos(alpha), -sin(30 + alpha), -sin(30 - alpha), -sin(30 - alpha), -sin(30 + alpha)])
    F_y_m = sin(theta) * np.array(6 * [-1]) #Creates an array of 6 * the item -cos(theta)
    #Moment around the x and z-axis (right hand rule). Since moments are around the point in the plane of the force application, only the vertical component matters
    M_x_m = sin(theta) * np.array([r, r, -0.5 * r, -0.5 * r, -0.5 * r, -0.5 * r])
    M_z_m = sin(theta) * np.array([0, 0, -cos(30) * r, -cos(30) * r, cos(30) * r, cos(30) * r])
    #Moment around the y-axis (downwards). Positive according to right hand rule
    M_y_m = cos(theta) * np.array(3 * [r * sin(alpha), -r * sin(alpha)])

    coeff_matrix = np.vstack([F_x_m, F_z_m, F_y_m, M_x_m, M_z_m, M_y_m])
    solution_matrix = np.array([F_x, F_z, F_y, M_x, M_z, M_y])
    forces = np.linalg.solve(coeff_matrix, solution_matrix)
    return(forces)
    
def Max_force(requirements, dimensions):
    #Set the max and min to the furthest infinities, to make sure any max or min applied later will ignore the initial value, without causing a crash or anything
    f_min = math.inf
    f_max = -math.inf
    #Set the initial non-varying parameters
    mass = dimensions["Mass"]
    safety_factor = requirements["Safety Factor"]
    #Test the forces for both the acceleration, and the deceleration case.
    for accel_y in (requirements["Longitudinal Acceleration"], requirements["Longitudinal Deceleration"]):
        requirements["F_y"] = accel_y * 9.81 * mass * safety_factor
        f_lat = requirements["Lateral Acceleration"] * 9.81 * mass * safety_factor
        #Due to symmetry (both rotational for every 120 degrees; and planar symmetry for every 60 degrees), the range of angles for forces to be tested is only  60 degrees. Any other angles would not give any new results, except that a different rod experiences the force, which doesn't change the final results.
        for angle in np.arange(0, 60, 0.1):
            requirements["F_x"] = f_lat * cos(angle)
            requirements["F_z"] = f_lat * sin(angle)
            forces = -Forces(requirements, dimensions) #Invert the forces, as in Forces, a positive force corresponds to compression inside the beam, but for stresses, a positive stress indicates tension
            #Find the new minimum and maximum forces.
            f_min = min(f_min, min(forces))
            f_max = max(f_max, max(forces))
    return(f_min, f_max)

if __name__ == "__main__":
    optimum = Find_optimal(materials, dimensions, requirements)
