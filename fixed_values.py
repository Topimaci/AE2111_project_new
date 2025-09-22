import numpy as np
#_______important for range calculations______
R_des = 6745000
f_cont = 0.05
e_f = 43000000
B = 4
v_cr = 200.62877
h_cr = 12500
t_E = 45*60
R_div = 277800


#______important for drag calculations_________
S_wet_over_S_w = 5.85
C_f = 0.004
psi = 0.0075
phi = 0.97
AR = 10


#________important for mass calculations_____
m_oe_over_m_MTO_fraction = 0.6077
m_pl_design = 750
m_pl_maxaximum = 1010


#______CRuise speed calculations_____
thrust_lapse = 0.24
beta = 0.94

#________Matching diagram: Landing and take-off distance_____
mass_fraction_landing = 0.7
landing_field = 950
density_landing = 1.165
C_L_max_landing = 2.3
density_takeoff = 1.225225
CL_takeoff = 1.6370864382454
temp_takeoff = 288.18
takeoff_field = 1250
oswald_efficiency = 0.887

#________Wing Loading______

wing_loading = np.arange(0,9100,100)


#________Matching diagram: climb rate______
C_d0 = 0.059
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

C_d0_119 = 0.0707175
C_d0_121a = 0.0447175
C_d0_121b = 0.0347175
C_d0_121c = 0.0087175
C_d0_121d = 0.0607175

e_119 = 0.901657174
e_121a = 0.849657174
e_121b = 0.849657174
e_121c = 0.797657174
e_121d = 0.901657174

