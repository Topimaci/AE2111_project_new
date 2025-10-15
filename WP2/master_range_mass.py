from Range_calculations import *
from mass_calculations import *
from dynamic_variables import *



###________________Calling First Range funcitons________________####
R_lost = R_lost_function(L_over_D_max, h_cr, v_cr)
R_eq_res = R_eq_res_function(R_div, t_E, v_cr)

R_eq = R_eq_function(R_des, R_lost, f_cont, R_eq_res)
R_aux = R_aux_function(R_des, R_eq)


#____________________Calling EXTRA VARIABLE FUNCTIONS_____________
TSFC = TSFC_function(B)
eta_j = eta_j_function(e_f, v_cr, TSFC)

fuel_mass_fraction = fuel_mass_fraction_function(R_eq, eta_j, e_f, L_over_D_max) #gets imported into mass_calculations

#____________________________________________________________#

m_MTO = m_MTO_function(m_pl_design, fuel_mass_fraction, m_oe_over_m_MTO_fraction)
m_oe = m_oe_function(m_oe_over_m_MTO_fraction, m_MTO)
m_f_des = m_f_des_function(m_MTO, m_pl_design, m_oe)
m_f_ferry = m_f_ferry_function(m_MTO, m_oe)                    #gets imported into range_calculations
m_f_harmonic = m_f_harmonic_function(m_MTO, m_oe, m_pl_maxaximum) #gets imported into range_calculations

print(m_MTO)

R_harmonic = R_harmonic_function(eta_j, L_over_D_max, e_f, m_f_harmonic, m_MTO, R_aux)
R_ferry = R_ferry_function(eta_j, L_over_D_max, e_f, m_f_ferry, m_MTO, R_aux)

if R_harmonic >= 6100000:
    print("Harmonic Range: PASS")
else: print("Harmonic Range: FAIL")

if R_ferry >= 7000000:
    print("Ferry Range: PASS")
else: print("Ferry Range: FAIL")

W_to = m_MTO * 9.81
new_S_w = W_to/designws
print(new_S_w)