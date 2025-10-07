import numpy as np
import math as m

#_______important for range calculations______
R_des = 6745000
f_cont = 0.05
e_f = 43000000
B = 4
v_cr = 200.63
h_cr = 12500
t_E = 45*60
R_div = 277800


#______important for drag calculations_________
S_wet_over_S_w = 5.85
C_f = 0.004
S_w = 43.95
S_wet = S_w * 5.85
psi = 0.0075
phi = 0.97
AR = 10


#________important for mass calculations_____
m_oe_over_m_MTO_fraction = 0.6077
m_pl_design = 750
m_pl_maxaximum = 1010


#______CRuise speed calculations_____
thrust_lapse = 0.24
beta = 0.95

#________Matching diagram: Landing and take-off distance_____
mass_fraction_landing = 0.7
landing_field = 700
density_landing = 1.165
C_L_max_landing = 1.51
density_takeoff = 1.225225
CL_takeoff = 1.6370864382454
temp_takeoff = 288.18
takeoff_field = 1250
oswald_efficiency = 0.887

#________Wing Loading______

wing_loading = np.arange(100,9100,100)
wing_loading_cs = np.arange(1200, 9100, 100)


#________Matching diagram: climb rate______
C_d0 = S_wet_over_S_w * C_f
p_ISO = 101325
T_ISO = 273.15
rho_ISO = 1.225225
## needs to be specified:
mass_fraction_climb = 0.97
climb_rate_requirement = 15  # m/s

#________Matching diagram: climb gradients______
cg_119 = 3.2
cg_121a = 0
cg_121b = 2.4
cg_121c = 1.2
cg_121d = 2.1


flap_to = 12
flap_la = 30


C_d0_119 = C_d0 + 0.02 +flap_la*0.0013
C_d0_121a = C_d0 + 0.02 +flap_to*0.0013
C_d0_121b = C_d0 + flap_to*0.0013
C_d0_121c = C_d0
C_d0_121d = C_d0 + flap_la*0.0013

e_119 = oswald_efficiency + flap_la * 0.0026
e_121a = oswald_efficiency + flap_to * 0.0026
e_121b = oswald_efficiency + flap_to * 0.0026
e_121c = oswald_efficiency
e_121d = oswald_efficiency + flap_la * 0.0026

mass_fraction_119 = 1
mass_fraction_121a = 1
mass_fraction_121b = 1
mass_fraction_121c = 1
mass_fraction_121d = 0.8

C_l_at_max_climb_gradient_119 = m.sqrt(C_d0_119 * m.pi * AR * e_119)
C_l_at_max_climb_gradient_121a = m.sqrt(C_d0_121a * m.pi * AR * e_121a)
C_l_at_max_climb_gradient_121b = m.sqrt(C_d0_121b * m.pi * AR * e_121b)
C_l_at_max_climb_gradient_121c = m.sqrt(C_d0_121c * m.pi * AR * e_121c)
C_l_at_max_climb_gradient_121d = m.sqrt(C_d0_121d * m.pi * AR * e_121d)

C_d_at_max_climb_gradient_119 = C_d0_119+C_l_at_max_climb_gradient_119**2/(e_119*AR*m.pi)
C_d_at_max_climb_gradient_121a = C_d0_121a+C_l_at_max_climb_gradient_121a**2/(e_121a*AR*m.pi)
C_d_at_max_climb_gradient_121b = C_d0_121b+C_l_at_max_climb_gradient_121b**2/(e_121b*AR*m.pi)
C_d_at_max_climb_gradient_121c = C_d0_121c+C_l_at_max_climb_gradient_121c**2/(e_121c*AR*m.pi)
C_d_at_max_climb_gradient_121d = C_d0_121d+C_l_at_max_climb_gradient_121d**2/(e_121d*AR*m.pi)