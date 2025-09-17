from Planform_DESIGN1 import *
from Drag_calculations import *
import Matching_Diagram.Minimum_speed as ms
import Matching_Diagram.Climb_rate as cr
import Matching_Diagram.Cruise_speed as cs
import Matching_Diagram.Landing_field_length as lfl
import Matching_Diagram.climb_grad as cg
import Matching_Diagram.Take_off_distance as td
import fixed_values as fv

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


#_________Wing Loading Calculations______________________________________________________________________________________________________
loads_minimum_speed = ms.Minimum_speed_function(fv.wing_loading)
loads_landing_field_length = lfl.landing_field_length(fv.mass_fraction, 950, 1.225, tudjaafaszom)
loads_cruise_speed = cs.cruise_speed(fv.beta,fv.thrust_lapse,fv.wing_loading,nemtom, 1.225, fv.v_cr, fv.AR, fv.e_f)
loads_climb_rate = 
loads_climb_grad_119 = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_119, fv.cd_119, fv.e_119, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121a = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_119, fv.cd_119, fv.e_119, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121b = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_119, fv.cd_119, fv.e_119, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121c = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_119, fv.cd_119, fv.e_119, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121d = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_119, fv.cd_119, fv.e_119, fv.AR, 1.225, nemtom, fv.B)
loads_to_field = td.take_off_distance()

print("Hello", sweep, taper, b , c_r, c_t, c_MAC, dihedral)
print(c_d0)
print(c_d_0_initial)
print(e_initial)
print(L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D )



