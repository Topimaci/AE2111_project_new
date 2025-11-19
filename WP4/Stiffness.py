import math as m
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

b = 19.585    # hard coded for now, should probably be pulled from somewhere in the code later on
max_displ = 0.15 * b
max_tip_rotat_deg = 10   # in degrees
max_tip_rotat_rad = m.radians(max_tip_rotat_deg)      # in radians
E = 71 * 10 ** 9    # Young's modulus
G = 27 * 10 ** 9    # Shear modulus

M_x =   # Import Moment function of M(y)
T =     # Import torque distribution function

y_breaks = #stringer breaks as np.array([...])
stringer_top_num = #nummber of stringer on these intervals np.array([...])
stringer_bottom_num = #same thing

#Linear interpolation of the stringers
string_top_interp = interp1d(y_breaks, stringer_top_num, kind="linear",
    fill_value="extrapolate")
string_bottom_interp = interp1d(y_breaks, stringer_bottom_num, kind="linear",
    fill_value="extrapolate")

y = sp.symbols("y")

spar_list = [lambda y: 0.4 * y + 0.1, 0.5, 1] # functions should be replaced, this is just an example, 0.5 is how much of the wing span the spar takes, 1 is how much of the chord it takes, measured from left side


def stiffness_distribution(y_pos, h_fs, h_rs, c_upper, c_lower, t, A_string, spar_list):
    # I Moment of Inertia Calculations
    #neutral axis
    x_c = (h_rs ** 2 + h_fs ** 2 + h_fs * h_rs) / (3 * (h_rs + h_fs))
    #Spar inertias
    I_fs = 1/12 * h_fs ** 3 * t + h_fs * t * (x_c - h_fs/2)**2
    I_rs = 1/12 * h_rs ** 3 * t + h_rs * t * (x_c - h_rs/2)**2
    #Skin inertias
    I_top = t * c_upper * (t/2 - x_c)**2
    I_bottom = t * c_lower * (((h_fs - x_c) + (h_rs - x_c))/2) ** 2

    num_top = float(string_top_interp(y_pos))
    num_bottom = float(string_bottom_interp(y_pos))

    #stringer inertias
    I_string_top = (A_string * (t - x_c)**2) * num_top
    I_string_bottom = (A_string * (((h_fs - x_c) + (h_rs - x_c))/2) ** 2) * num_bottom
       
    if spar_list != []:
        I_step = 0
        for h_spar_func, y_crit in spar_list:

            # h_spar is a sympy expression of y
            h_spar_y = h_spar_func(y)
            # compute I_spar(y)
            I_spar_y = 1/12 * h_spar_y**3 * t
            # add step contribution
            I_step += sp.Piecewise(
                (I_spar_y, y < y_crit),
                (0, True)
            )
        I_total = I_step + I_string_bottom + I_string_top + I_bottom + I_top + I_fs + I_rs
        a = spar_list[2]
        w = c_upper - a 
        lefthand_matrix = np.array([[(2*w+2*a), -w, -2*a*w*G*t], [-w, 4*b, - 2*w**2*G*t], [2*a*w, 2 * w**2, 0]])
        righthand_matrix = np.array([0, 0, 1])
        solution = np.linalg.solve(lefthand_matrix, righthand_matrix)
        q1, q2, dtheta_dy = solution
        J = 1 / (G * dtheta_dy)

    else:
        I_total = I_string_bottom + I_string_top + I_bottom + I_top + I_fs + I_rs
        A = (h_fs + h_rs) / 2 * c_upper
        circ = 1/ t * (h_fs + c_upper + h_rs + c_lower)
        J = 4 * A ** 2 / circ
    return I_total, J

I_xx, J = stiffness_distribution()
d2v_dy2 = - M_x / (E * I_xx)
dth_dy = T / (G * J)

estimate_dv, error_dv = sp.integrate.quad(d2v_dy2, 0, b)
estimate_v, error_v = sp.integrate.quad(estimate_dv, 0, b)
estimate_th, error_th = sp.integrate.quad(dth_dy, 0, b)

