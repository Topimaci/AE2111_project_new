import math as m


def calculate_geometric_parameters_wing(S_w, AR, M_cr):

    sweep = m.acos(1.16/ (M_cr + 0.5) ) * 180 / m.pi

    taper = 0.2 *(2 - sweep * m.pi / 180 )

    b = m.sqrt(AR*S_w)

    c_root = (2 * S_w)/((1+taper) * b)

    c_tip = c_root * taper

    c_MAC = (2/3)* c_root * ((1 + taper + taper**2)/(1+taper))
    
    dihedral = 3 - 0.1 * sweep + 2 #3 default for unswept, subtract 0.1 for every degree of sweep, +2 from low wing config

    sweep_LE = m.atan((-1/4*c_tip + 1/4*c_root + m.tan(sweep)*b) / b ) *180 /m.pi

    return sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE


def calculate_aerodynamic_performance(thickness_to_chord):

    c_d0 = 0.0035 + 0.018 * thickness_to_chord

    return c_d0


def calculate_MAC_position(b, c_root, c_tip, sweep):

    b = b/2 #we used half the wingspan because in the geometric 
    print(b)

    sweep = sweep * m.pi/180
    print(c_root, c_tip, sweep)


    a_1 = (b/(5/4 * c_tip + 7/4 * c_root - m.tan(sweep)*b))**(-1)

    a_2 = (b/(3/4 *(c_root + c_tip) - m.tan(sweep) * b - c_tip - 0.5 * c_root ) )**(-1)

    a_3 = (5/4 * c_tip + 3/4 * c_root - m.tan(sweep)* b - c_tip - c_root)/b

    a_4 = (1/4 * c_tip + 3/4 * c_root - m.tan(sweep)* b - c_tip)/b

    b_2 = c_tip + 0.5 * c_root

    b_3 = c_tip + c_root

    b_4 = c_tip

    y_spanwise = b_2 / (a_1 - a_2)

    xlemac = -a_3 * y_spanwise

    lengthMAC = -xlemac + b_3 - a_4 * y_spanwise - b_4



    return y_spanwise, xlemac, lengthMAC




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
# - dihedral = dihedral angle (in degrees)

# Assumptions:
# - MGC roughly MAC
# - low wing config, for dihedral


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
