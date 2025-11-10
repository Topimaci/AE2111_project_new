import math as m
import io
import contextlib
import os
import sys

# Capture any print output during imports
f = io.StringIO()
with contextlib.redirect_stdout(f):
    # Current folder (functions/)
    current_dir = os.path.dirname(__file__)
    
    # Add the project root (parent of functions, WP2, variables) to sys.path
    project_root = os.path.abspath(os.path.join(current_dir, ".."))
    sys.path.insert(0, project_root)

    # Now you can use package-style imports
    import WP2.HLDs as HLDs
    import WP2.master_range_mass as mrm
    import variables.fixed_values as fv
    import functions.Planform_DESIGN1 as Pl


## Required values for roll performance
c_l_alpha = 0.11965 *360 / 2 * m.pi ## per radian
S_w = 30.48
taper, b, c_root, c_tip, c_MAC, dihedral = Pl.calculate_geometric_parameters_wing(S_w, fv.AR, 0.68)
k_const = 2.5 ## misc constant
m_const = 1.8 ## misc constant
rear_spar_pos = 0.65 ## as fraction of chord
c_d0 = Pl.calculate_aerodynamic_performance(0.12)
deflection_up = 20 ## degrees
deflection_down = 0.75 * deflection_up ## degrees
deflection = 1/2 * (deflection_up + deflection_down) ## degrees
gap = 0.02 ## fraction of semispan
V = 1.23 * m.sqrt(2 * mrm.m_MTO * 9.80665 / (fv.rho_ISO * S_w * 1.37))
roll_performance_requirement = 60/11 ## degrees per second

b_1 =  (HLDs.end_pos_span_TE + gap) * b/2  ## in meters
db = 0.00001


## Roll damping coefficient
def C_lP(c_l_alpha, c_d0, S_w, b, c_tip, c_root):
    return - (c_l_alpha + c_d0)/S_w * b/8 * (c_tip + 1/3 * c_root)

## Aileron control derivative
def C_ldalpha(c_l_alpha, S_w, b, k_const, m_const, rear_spar_pos, c_tip, c_root, b_1, b_2):
    return 2 * c_l_alpha / (S_w * b) * k_const * (1 - rear_spar_pos) / (1 + m_const * (1 - rear_spar_pos)) * \
        ((2/b * (c_tip - c_root) * b_2**3 / 3 + c_root * b_2**2 / 2) - \
         (2/b * (c_tip - c_root) * b_1**3 / 3 + c_root * b_1**2 / 2))

## Rate of roll
def roll_performance(deflection, V, b, b_2): ## degrees per second
    return - C_ldalpha(c_l_alpha, S_w, b, k_const, m_const, rear_spar_pos, c_tip, c_root, b_1, b_2) /\
        C_lP(c_l_alpha, c_d0, S_w, b, c_tip, c_root) * deflection * 2 * V / b

## Lengthening the aileron until the roll performance requirement is met
def b_2_repetition(b_1, db):
    smaller = True
    b_2 = b_1
    while smaller:
        b_2 = b_2 + db
        roll_performance_calculated = roll_performance(deflection, V, b, b_2)
        if roll_performance_calculated > roll_performance_requirement:
            smaller = False
            return b_2
        if b_2 > b/2:
            smaller = False
            return 0


## print(b_2_repetition(b_1, db)/(b/2))  ## as fraction of semi-span
## print(roll_performance(deflection, V, b, 0.95*(b/2)))  ## degrees per second