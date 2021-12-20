import math

def Euler_column(material, dimensions, thickness):
    """
    :param material: Dict containing material properties
    :param dimensions: Dict containing the dimensions of the tank
    
    Returns the critical euler column buckling stress, derived from the given dimensions and material properties.
    WARNING: Does NOT yet include any safety factors.
    """
    #Extract the parameters from the dicts
    length = dimensions["Length"]
    radius = dimensions["Radius"]

    E = material["Modulus"]
    #Calculate the material properties.
    I = 1 / 4 * math.pi * ((radius + thickness/2)**4 - (radius - thickness/2)**4)
    A_cross = Area(radius, thickness) #Calculate the cross-sectional area.
    sigma_crit = (math.pi * math.pi * E * I) / (A_cross * length * length)
    return(sigma_crit)


def Shell_buckling(material, dimensions, thickness):
    """
    :param material: Dict containing material properties
    :param dimensions: Dict containing the dimensions of the tank
    :param requirements: Dict containing the requirements for the tank wall
    
    Returns the critical shell buckling stress, derived from the given dimensions and material properties".
    WARNING: Does NOT yet include any safety factors.
    """
    length = dimensions["Length"]
    radius = dimensions["Radius"]

    E = material["Modulus"]
    poisson = material["Poisson"]

    lambd = length**2 * math.sqrt(12 * (1 - poisson**2)) / (math.pi**2 * radius * thickness)
    k = 2 * lambd #Love this. Minimum value of y = lambd + a / lambd, where a == lambd^2

    sigma_crit = 1.983 * k * math.pi ** 2 * E / (12 * (1 - poisson**2)) * (thickness / length)**2
    return(sigma_crit)


def sort_by(iterator, index):
    lst = [[item[index], item] for item in iterator] #Create a new list, where each item has the value at the desired index as the first value
    lst = list(sorted(lst)) #Sort the list using sorted. Since the given index is the first value in each item, the list is sorted using this value
    lst = [item[1] for item in lst] #Remove the sorting criteria again, so the list is just a list of normal items again
    return(lst)

def Area(radius, thickness):
    """
    Returns the area of a cylinder with a given median radius and thickness.
    """
    return(math.pi * ((radius + thickness/2)**2 - (radius - thickness/2)**2))
