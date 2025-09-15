import math as m




def calculate_geometric_parameters_wing(S_w, AR, M_cr):

    sweep = m.acos(1.16/ (M_cr + 0.5) ) * 180 / m.pi

    taper = 0.2 *(2 - sweep * m.pi / 180 )

    b = m.sqrt(AR*S_w)

    c_root = (2 * S_w)/((1+taper) * b)

    c_tip = c_root * taper

    c_MAC = (2/3)* c_root * ((1 + taper + taper^2)/(1+taper))

    return sweep, taper, b, c_root, c_tip, c_MAC


def calculate_aerodynamic_performance_based_on_given_geometric_parameters(thickness_to_chord, C_Lmax,):

    c_d0 = 0.0035 + 0.018 * thickness_to_chord

    c_lmax = 


    return c_d0, c_lmax


#----------------------------------------------------------------------------------------------------#
# FUNCTION: calculate_geometric_parameters_wing

# Inputs:

# S_w = Surface area wing
# AR = Aspect ratio
# M_cr = Mach number at cruise

# Outputs:

# - sweep = Sweep angle (in degrees)
#  - b = winsgpan
# - taper = taper ratio
# - c_root = root chord
# - c_tip = tip chord

# Assumptions:
# - MGC roughly MAC


#----------------------------------------------------------------------------------------------------#

# FUNCTION: calculate_geometric_parameters_wing

# Inputs:

# C_Lmax

# Outputs "calculate_aerodynamic_performance_based_on_given_geometric_parameters": 

# c_d0 
# c_lmax

# Assumptions "calculate_aerodynamic_performance_based_on_given_geometric_parameters":

# -x/c front spar assumed to be 0.2 (adsee book)
# -x/c rear spar assumed to be 0.7 (adsee book)
# thickness_to_chord is assumed to be within 0.06 and 0.25 for the formula from adsee to be valid 

#----------------------------------------------------------------------------------------------------#
