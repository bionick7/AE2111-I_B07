import math
from numpy import arange
import matplotlib.pyplot as plt

from math import exp, pi

materials = {"Aluminium-2024-T4": #https://www.gabrian.com/2024-aluminum-properties/
             {"Modulus": 73.1e9,
              "Yield strength": 324e6, #Pa
              "Ultimate strength": 469e6, #Pa
              "Density": 2780, #kg/m^3
                 },
             "Aluminium-7075-T6": #https://www.gabrian.com/7075-aluminum-properties/
             {"Modulus": 71.7e9,
              "Yield strength": 503e6,
              "Ultimate strength": 572e6,
              "Density": 2810,
                 },
             "Ti-6AI-4V": #http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MTP641
             {"Modulus": 113.8e9,
              "Yield strength": 880e6,
              "Ultimate strength": 950e6,
              "Density": 4430,
              }
             }


SC_parameters = {"Mass": 1903, #kg
                 "S/C length": 8.790, #m (from bottom of cylinder to S/C top
                 "Cylinder length": 2.5, #m
                 "CoM offset": 0, #kg*m; Should be re-calculated in each loop as (1/2 * diameter + 0.5) * 200
                 "Max cylinder diameter": 2.3, #m
                 }

requirements = {"Freq longitudinal": 31, #Hz
                "Freq lateral": 9, #Hz
                "Axial acceleration": 6, #G's
                "Axial deceleration": -2.5, #G's
                "Lateral acceleration": 2, #G's
                "Design safety factor": 1.25,
                "Yield safety": 1.1, #On top of the Design safety factor
                "Ultimate safety": 1.25, #On top of the Design safety factor
                }


def thickness_structure(material, SC_parameters, requirements, diameter) -> float:
    """
    Calculates the thickness required for the structural component of the S/C to
    withstand all loads (with safety factors). For frequency calculations, the
    program uses the "simple case", where all mass is at the end/top of the S/C.
    """
    #Unpack all variables
    E = material["Modulus"]
    density = material["Density"]
    sigma_yield = material["Yield strength"]
    sigma_ult = material["Ultimate strength"]

    SC_length = SC_parameters["S/C length"]
    cylinder_length = SC_parameters["Cylinder length"]
    mass = SC_parameters["Mass"]

    safety_factor = requirements["Design safety factor"]
    safety_yield = requirements["Yield safety"]
    safety_ult = requirements["Ultimate safety"]
    #Set all design loads ( / accelerations):
    CoM_offset = SC_parameters["CoM offset"] * safety_factor
    posG = requirements["Axial acceleration"] * safety_factor
    negG = requirements["Axial deceleration"] * safety_factor
    sideG = requirements["Lateral acceleration"] * safety_factor

    freq_long = requirements["Freq longitudinal"] * safety_factor
    freq_lat = requirements["Freq lateral"] * safety_factor

    #Calculate the design bending moments that occur during positive and
    #negative axial G loading and maximum lateral loading. If CoM_offset = 0,
    #they will be the same. Calculated values are all positive.
    M_posG = (mass * sideG * SC_length * 9.81) + (CoM_offset * posG * 9.81)
    M_negG = (mass * sideG * SC_length * 9.81) + (CoM_offset * abs(negG) * 9.81)
    
    #Calculate the longitudinal design loads for positive and negative axial
    #loading
    F_long_posG = posG * mass * 9.81
    F_long_negG = negG * mass * 9.81

    #Calculate the Sigmat's. Explanation of what the Sigmat are in this program, can be
    #found in the function "Buckling".
    Sigmat_posG_comp = (F_long_posG/(math.pi * diameter)) + (4 * M_posG / (math.pi * (diameter**2)))
    Sigmat_negG_comp = (F_long_negG/(math.pi * diameter)) + (4 * M_negG / (math.pi * (diameter**2)))
    Sigmat_posG_tens = (F_long_posG/(math.pi * diameter)) - (4 * M_posG / (math.pi * (diameter**2)))
    Sigmat_negG_tens = (F_long_negG/(math.pi * diameter)) - (4 * M_negG / (math.pi * (diameter**2)))

    t_buckling = buckling_structure(E, SC_length, diameter, safety_ult * max(Sigmat_posG_comp, Sigmat_negG_comp))

    #Find the value which will be the most critical for the wall thickness
    #Absolute values are taken, since yielding/failure is assumed to happen
    #at the same stress in both modes.
    Sigmat_max = max((abs(Sigmat_posG_comp), abs(Sigmat_negG_comp), abs(Sigmat_posG_tens), abs(Sigmat_negG_tens)))

    t_yield = Sigmat_max * safety_yield / sigma_yield
    t_ultimate = Sigmat_max * safety_ult / sigma_ult

    #Calculate the A*E required for the longitudinal frequency requirement
    AE_req_long = ((freq_long / 0.160)**2) * mass * SC_length
    #Calculate the E*I required for the lateral frequency requirement
    EI_req_lat = ((freq_lat / 0.276)**2) * mass * (SC_length**3)
    t_freq_long = AE_req_long / (math.pi * diameter * E)
    t_freq_lat = (8 * EI_req_lat) / (math.pi * (diameter**3) * E)
    return(max(t_buckling, t_yield, t_ultimate, t_freq_long, t_freq_lat))
                 

def buckling_structure(E, length, diameter, Sigmat_max_comp) -> float:
    """
    Iterates over thicknesses, until a thickness is found for which buckling
    will not occur. The iterative approach is used since separating thickness is
    quite difficult.
    """
    #Sigmat_max_comp is the maximum "Force" that will occur. While in this state, it
    #isn't a conventional unit, it can however quickly be turned into stress by
    #simply dividing by the cylinder thickness, and it is thus basicallly
    #stress * thickness
    thickness = 0
    thickness_step = 1e-6
    while True:
        thickness += thickness_step
        stress = Sigmat_max_comp / thickness
        critical_stress = (9 * (2 * thickness / diameter)**(1.6) + 0.16 * (thickness / length)**(1.3)) * E
        if stress <= critical_stress:
            return thickness


def mass_structure(material, diameter, thickness) -> float:
    """
    Returns the mass of the cylinder for the given configuration.
    """
    return(math.pi * diameter * thickness * material["Density"])


def structure_design_routine():
    step_size = 0.001
    max_diameter = SC_parameters["Max cylinder diameter"] + step_size #Set the diameter up to which the arange should go
    optimums = [] #Initialize the list for all optimum values. This list is usefull if anyone wants to automate the finding of the optimal material.
    for material_name in materials: #For each material:
        material = materials["Aluminium-2024-T4"]
        data = [] #Initialize the list of data points
        for diameter in arange(0.5, max_diameter, step_size):
            SC_parameters["CoM offset"] = ((diameter / 2) + 0.5) * 190 #Calculate the CoM offset of the orbiter
            thickness = thickness_structure(material, SC_parameters, requirements, diameter) #Get the required thickness for the given diameter
            data.append([diameter, mass_structure(material, diameter, thickness), thickness]) #Append the data point to the data list
        #Find the optimum values
        masses = [point[1] for point in data]
        M_min = min(masses)
        t_opt = data[masses.index(M_min)][2]
        d_opt = data[masses.index(M_min)][0]
        #Print the optimum values
        print(f"Material: {material_name}")
        print(f"Minimum mass achieved at a diameter of {d_opt:.3f} m, with thickness {t_opt:.5f} m. Mass is {M_min:.3f} kg.")
        #Add the optimum values to the optimums list
        optimums.append([material_name, M_min, t_opt, d_opt])
        #Add the diameter-mass plot to the plot
        plt.plot([point[0] for point in data], masses, label = material_name)
    #Finalise the plot
    plt.legend()
    plt.title("Mass vs cylinder diameter")
    plt.xlabel("Cylinder diameter [m]")
    plt.ylabel("Cylinder mass [kg]")
    plt.show()


def adcs_sizing(dry_mass):
    """
    dry_mass -> component_dry_mass, wet_mass
    """
    dv = 1442
    w = 2864
    m_prop = dry_mass * (exp(dv/w) - 1)
    m_tank = adcs_tank_mass(m_prop)
    return m_tank, m_prop + m_tank + dry_mass


def adcs_tank_mass(propellant_mass):
    return 0.0348 * propellant_mass + 58.15 - 10.92


def prop_sizing(wet_mass):
    """
    wet_mass -> component_dry_mass, launch mass
    """
    m_prop = 1000
    r = 1
    m_tank = 0
    while abs(r) > 0.001:
        if r > 0:
            m_prop *= 0.9
        else:
            m_prop *= 1.1
        r = prop_eq(m_prop, wet_mass - m_tank)
        m_tank = 0.1171 * m_prop ** 0.848
    return m_tank, m_prop + wet_mass


def prop_eq(x, pl):
    """
    needs to be zero. payload mass pl is a parameter
    """
    return (0.1171 * x**0.848 + x + pl) / (0.1171 * x**0.848 + pl) - exp(10_500/10_000)


def structure_sizing(launch_mass):
    """
    launch_mass -> component_dry_mass
    """
    sc_params = SC_parameters.copy()
    sc_params["Mass"] = launch_mass
    material = materials["Aluminium-2024-T4"]
    thickness = thickness_structure(material, sc_params, requirements, 2.3)
    print("thickness = " + str(thickness) + "m")
    return mass_structure(material, 2.3, thickness)


def convergence_check():
    dry_mass = 0
    wet_mass = 0
    launch_mass = 0
    constant_dry = 886.12  # kg From Budget before reverting comms
    adcs_dry, prop_dry, struc_dry = 0, 0, 0
    for i in range(100):
        dry_mass = constant_dry + adcs_dry + prop_dry + struc_dry
        adcs_dry, wet_mass = adcs_sizing(dry_mass)
        prop_dry, launch_mass = prop_sizing(wet_mass)
        struc_dry = structure_sizing(launch_mass)
    print(f"""
        Dry:    {dry_mass:>7.3f} kg
        Wet:    {wet_mass:>7.3f} kg
        Launch: {launch_mass:>7.3f} kg
        
        Dry mass Per component:
        ADCS Tanks {adcs_dry:>7.3f} kg
        H2 Tanks   {prop_dry:>7.3f} kg
        Structure  {struc_dry:>7.3f} kg
        
        Fuel masses:
        ADCS {wet_mass - dry_mass:>7.3f} kg
        H2   {launch_mass - wet_mass:>7.3f} kg
        H2 volume {(3/4/pi * (launch_mass - wet_mass)/70)**.3333}
        """)


if __name__ == '__main__':
    convergence_check()
