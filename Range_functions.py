import math


def fuel_mass_fraction_function(R_eq, eta_j, e_f, L_over_D_max):
    fuel_mass_fraction = 1 - math.exp(-R_eq/(eta_j*(e_f/9.81)*L_over_D_max))
    return fuel_mass_fraction

def R_lost_function(L_over_D_max, h_cr, v_cr):
    R_lost = 1/0.7 *L_over_D_max*(h_cr+(v_cr**2)/(2*9.81))
    return  R_lost

def R_eq_function(R_des, R_lost, f_cont, R_eq_res):
    R_eq = (R_des + R_lost)*(1+f_cont) + R_eq_res
    return R_eq

def eta_j_function(e_f, v_cr, TSFC):
    eta_j = v_cr/(TSFC*e_f)
    return eta_j

def 