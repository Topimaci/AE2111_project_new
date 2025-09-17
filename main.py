from Planform_DESIGN1 import *
from Drag_calculations import *

thickness_to_chord = 0.15 # assumption 
AR = 10
M_cr = 0.68   #this is for requirement, maybe 0.7 for cruise?
S_w = 30.46

#for drag calculations the values below are assumed
S_wet_over_S_w = 5.85
C_f = 0.004
psi = 0.0075
phi = 0.97


sweep, taper, b , c_r, c_t, c_MAC, dihedral = calculate_geometric_parameters_wing(S_w, AR, M_cr)
c_d0 = calculate_aerodynamic_performance(thickness_to_chord)

c_d_0_initial = c_d0_function(S_wet_over_S_w, C_f)

e_initial = e_function(psi, phi, AR)

L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D = L_over_D_max_function(AR,e_initial,c_d_0_initial)





print("Hello", sweep, taper, b , c_r, c_t, c_MAC, dihedral)
print(c_d0)
print(c_d_0_initial)
print(e_initial)
print(L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D )



