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

#________Wing Loading______

wing_loading = np.arange(0,9100,100)
