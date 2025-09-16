import math

from Drag_calculations import L_over_D_max


R_des = 6745000
f_cont = 0.05
e_f = 43000000
B = 4
v_cr = 200.62877
h_cr = 12500
t_E = 45*60
R_div = 277800

##_____________first range functions___________________________--##
def R_lost_function(L_over_D_max, h_cr, v_cr):
    R_lost = 1/0.7 *L_over_D_max*(h_cr+(v_cr**2)/(2*9.81))
    return  R_lost

def R_eq_res_function(R_div, t_E, v_cr):
    R_eq_res = 1.2*R_div + t_E*v_cr
    return R_eq_res

def R_eq_function(R_des, R_lost, f_cont, R_eq_res):
    R_eq = (R_des + R_lost)*(1+f_cont) + R_eq_res
    return R_eq

def R_aux_function(R_des, R_eq):
    R_aux = R_eq - R_des
    return  R_aux
##____________________________________________________________________________###


#______________________EXTRA VARIABLE CALCULATIONS________________________
def TSFC_function(B):
    TSFC = 22*(B**(-0.19))
    return TSFC

def eta_j_function(e_f, v_cr, TSFC):
    eta_j = v_cr/(TSFC*e_f)
    return eta_j
#_______________________________________________________________



#________________SECOND RANGE CALCULATIONS______________________________
from mass_calculations import m_f_ferry
from mass_calculations import m_f_harmonic

def fuel_mass_fraction_function(R_eq, eta_j, e_f, L_over_D_max): ####this functions gets imported into mass_calculations
    fuel_mass_fraction = 1 - math.exp(-R_eq/(eta_j*(e_f/9.81)*L_over_D_max))
    return fuel_mass_fraction


def R_harmonic_function(eta_j, L_over_D_max, e_f, m_f_harmonic, m_MTO, R_aux):
    R_harmonic = eta_j * L_over_D_max * (e_f/9.81)*math.log(1-m_f_harmonic/m_MTO) - R_aux
    return R_harmonic

def R_ferry_function(eta_j, L_over_D_max, e_f, m_f_ferry, m_MTO, R_aux):
    R_ferry = eta_j * L_over_D_max * (e_f/9.81)*math.log(1-m_f_ferry/m_MTO) - R_aux
    return R_ferry
#__________________________________________________________________________-



###________________Calling First Range funcitons________________####
R_lost = R_lost_function(L_over_D_max, h_cr, v_cr)
R_eq_res = R_eq_res_function(R_div, t_E, v_cr)

R_eq = R_eq_function(R_des, R_lost, f_cont, R_eq_res)
R_aux = R_aux_function(R_des, R_eq)
#____________________________________________________________#


#____________________Calling EXTRA VARIABLE FUNCTIONS_____________
TSFC = TSFC_function(B)
eta_j = eta_j_function(e_f, v_cr, TSFC)

#______________________________________________________________-

#_____________CALLING SECOND RANGE FUNCTIONS____________________
fuel_mass_fraction = fuel_mass_fraction_function(R_eq, eta_j, e_f, L_over_D_max) #gets imported into mass_calculations
R_harmonic = R_harmonic_function(eta_j, L_over_D_max, e_f, m_f_harmonic)
R_ferry = R_ferry_function(eta_j, L_over_D_max, e_f, m_f_ferry)

#_______________________________________

print(R_harmonic)