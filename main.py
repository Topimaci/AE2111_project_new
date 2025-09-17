from Planform_Design import *
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


print("Hello", sweep, taper, b , c_r, c_t, c_MAC, dihedral)
print(c_d0)



