import math as m
import fixed_values as fv
import Planform_DESIGN1 as Pl

c_l_alpha = 2 * m.pi  ## modifiable
S_w = fv.S_w
sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE = Pl.calculate_geometric_parameters_wing(S_w, fv.AR, 0.68, 24.5)
k = 2.5 ## misc constant
m = 1.8 ## misc constant
rear_spar_pos = 0.7 ## as fraction of chord
c_d0 = Pl.calculate_aerodynamic_performance(0.12)
deflection_up = 20 ## degrees
deflection_down = 0.75 * deflection_up ## degrees
deflection = 1/2 * (deflection_up + deflection_down)
V_cruise = fv.v_cr
roll_performance_requirement = 45/1.4 ## degrees per second for class II

b_1 = 0.85 * b/2 ## start of aileron
b_2 = 0.9 * b/2 ## end of aileron

def C_lP(c_l_alpha, c_d0, S_w, b, c_tip, c_root):
    return - (c_l_alpha + c_d0)/S_w * b/8 * (c_tip + 1/3 * c_root)

def C_ldalpha(c_l_alpha, S_w, b, k, m, rear_spar_pos, c_tip, c_root, b_1, b_2):
    return 2 * c_l_alpha / (S_w * b) * k * rear_spar_pos / (1 + m * rear_spar_pos) * \
        ((2/b * (c_tip - c_root) * b_2**3 / 3 + c_root * b_2**2 / 2) - \
         (2/b * (c_tip - c_root) * b_1**3 / 3 + c_root * b_1**2 / 2))

def roll_performance(deflection, V_cruise, b): ## degrees per second
    return - C_ldalpha(c_l_alpha, S_w, b, k, m, rear_spar_pos, c_tip, c_root, b_1, b_2) /\
        C_lP(c_l_alpha, c_d0, S_w, b, c_tip, c_root) * deflection * 2 * V_cruise / b

print(roll_performance(deflection, V_cruise, b))