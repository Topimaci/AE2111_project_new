import math as m
import sympy as sp

b = 19.585    # hard coded for now, should probably be pulled from somewhere in the code later on
max_displ = 0.15 * b
max_tip_rotat_deg = 10   # in degrees
max_tip_rotat_rad = m.radians(max_tip_rotat_deg)      # in radians

y = sp.symbols("y")

spar_list = [(0.4, 0.5)]


def stiffness_distribution(h_fs, h_rs, c_upper, c_lower, t, A_string, num_string_top, num_string_bottom, spar_list):
    x_c = (h_rs ** 2 + h_fs ** 2 + h_fs * h_rs) / (3 * (h_rs + h_fs))
    I_fs = 1/12 * h_fs ** 3 * t + h_fs * t * (x_c - h_fs/2)**2
    I_rs = 1/12 * h_rs ** 3 * t + h_rs * t * (x_c - h_rs/2)**2
    I_top = t * c_upper * (t/2 - x_c)**2
    I_bottom = t * c_lower * (((h_fs - x_c) + (h_rs - x_c))/2) ** 2
    I_string_top = (A_string * (t - x_c)**2) * num_string_top
    I_string_bottom = (A_string * (((h_fs - x_c) + (h_rs - x_c))/2) ** 2) * num_string_bottom
    
    Step = 0
    I_spar = 1/12 * spar_list[0] ** 3 * t
    for I_spar 

# Loop through all contributions and add as Piecewise terms
for I_extra, y_crit in step_list:
    I_step += sp.Piecewise(
        (I_extra, y < y_crit),
        (0, True)
    )
    I_spar = sp.Piecewise((I_spar1, y<0.5), (0, True))
    
    
