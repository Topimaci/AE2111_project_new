import math as m
import numpy as np


def calculate_geometric_parameters_wing(S_w, AR, M_cr):
    """
    Calculate Geometric Parameters of the Wing

    :param S_w: Wing Surface Area [m^2]
    :param AR: Aspect Ratio [-]
    :param M_cr: Mach Cruise Number [-]
    :return: Tuple (Sweep Angle for c/4 [deg], Taper Ratio, Wingspan [m], Root Chord [m], Tip Chord [m], Dihedral Angle [deg], Sweep Angle for LE [deg])
    """
    sweep = m.acos(1.16/ (M_cr + 0.5) ) * 180 / m.pi
    
    taper = 0.2 *(2 - sweep * m.pi / 180 )

    b = m.sqrt(AR*S_w)

    c_root = (2 * S_w)/((1+taper) * b)

    c_tip = c_root * taper

    c_MAC = (2/3)* c_root * ((1 + taper + taper**2)/(1+taper))
    
    dihedral = 3 - 0.1 * sweep + 2 #3 default for unswept, subtract 0.1 for every degree of sweep, +2 from low wing configuration

    sweep_LE = m.atan((-1/4*c_tip + 1/4*c_root + m.tan(sweep * m.pi / 180)*b) / b ) *180 /m.pi

    return  taper, b, c_root, c_tip, c_MAC, dihedral


def C_L_design(M_MTO, v_cruise, density_cruise, Wing_area):
    C_L_design = M_MTO*9.81/(0.5* v_cruise**2 * density_cruise * Wing_area)
    return  C_L_design


def sweep_drag_divergence(C_L):
    coeffs = [0.68, -0.87, 0.14, C_L]

    roots = np.roots(coeffs)

    # Find the real root where -1 <= cos(Λ) <= 1
    for r in roots:
        if abs(r.imag) < 1e-9:
            x = r.real
            if -1 <= x <= 1:
                return np.degrees(np.arccos(x))  # Λ in degrees

    return None  # no physical solution found






#def calculate_aerodynamic_performance(thickness_to_chord):   #######newer estimations for Cd0
    """
    Calculate Aerodynamic Performance Based on Thickness-To-Chord Ratio

    :param thickness_to_chord: Thickness-To-Chord Ratio [-]
    :return: Parasitic Drag Coefficient c_d0 [-]
    """
    #c_d0 = 0.0035 + 0.018 * thickness_to_chord

    #return c_d0


def calculate_MAC_position(b, c_root, c_tip, sweep):
    """
    Calculate MAC Position

    :param b: Wingspan [m]
    :param c_root: Root Chord [m]
    :param c_tip: Tip Chord [m]
    :param sweep: Sweep Angle for c/4 [deg]
    :return: Tuple (y-spanwise location of MAC [m], x-spanwise location of MAC [m], Chord length of MAC [m])
    """

    b = b/2 #For Calculations We Used b as half-wingspan
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



    return y_spanwise, xlemac

def sweep_converter(sweep_LE, chord_root, taper, x_over_c, span):
    sweep_x_over_c = m.atan(m.tan(sweep_LE*m.pi/180) - x_over_c * (2*chord_root/span)*(1-taper)) * 180/m.pi
    return sweep_x_over_c 

    

