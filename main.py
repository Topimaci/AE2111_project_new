from Planform_DESIGN1 import *
from Drag_calculations import *
import Matching_Diagram.Minimum_speed as ms
import Matching_Diagram.Climb_rate as cr
import Matching_Diagram.Cruise_speed as cs
import Matching_Diagram.Landing_field_length as lfl
import Matching_Diagram.climb_grad as cg
import Matching_Diagram.Take_off_distance as td
import fixed_values as fv
import matplotlib.pyplot as plt

thickness_to_chord = 0.15 # assumption 
AR = 10
M_cr = 0.68   #this is for requirement, maybe 0.7 for cruise?
S_w = 30.46

#for drag calculations the values below are assumed
S_wet_over_S_w = 5.85
C_f = 0.004
psi = 0.0075
phi = 0.97


sweep, taper, b , c_r, c_t, c_MAC, dihedral, sweep_LE = calculate_geometric_parameters_wing(S_w, AR, M_cr)

c_d0 = calculate_aerodynamic_performance(thickness_to_chord)

c_d_0_initial = c_d0_function(S_wet_over_S_w, C_f)

e_initial = e_function(psi, phi, AR)

L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D = L_over_D_max_function(AR,e_initial,c_d_0_initial)

y_spanwise, xlemac, lengthMAC = calculate_MAC_position(b, c_r, c_t, sweep)

#_________Wing Loading Calculations______________________________________________________________________________________________________
loads_minimum_speed = ms.Minimum_speed_function(fv.wing_loading)
loads_landing_field_length = lfl.landing_field_length(fv.mass_fraction, 950, 1.225, fv.C_L_max_landing)
loads_cruise_speed = cs.cruise_speed(fv.beta,fv.thrust_lapse,fv.wing_loading,fv.C_d0, 1.225, fv.v_cr, fv.AR, fv.e_f)
loads_climb_rate = cr.climb_rate(fv.wing_loading)
loads_climb_grad_119 = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_119, fv.C_d0_119, fv.e_119, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121a = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_121a,fv.C_d0_121a, fv.e_121a, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121b = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_121b, fv.C_d0_121b, fv.e_121b, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121c = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_121c, fv.C_d0_121c, fv.e_121c, fv.AR, 1.225, nemtom, fv.B)
loads_climb_grad_121d = cg.climb_grad(fv.wing_loading, massfraction, fv.cg_121d, fv.C_d0_121d, fv.e_121d, fv.AR, 1.225, nemtom, fv.B)
V2 = td.find_V_2(fv.wing_loading, fv.density_takeoff, fv.CL_takeoff)
Mach = td.find_Mach(V2, fv.temp_takeoff)
theta_delta = td.find_theta_delta(fv.temp_takeoff, Mach)
alpha = td.find_alpha_t(theta_delta, Mach, fv.B)
loads_to_field = td.take_off_distance(alpha, fv.wing_loading, fv.takeoff_field, fv.density_takeoff, fv.oswald_efficiency, fv.AR)

print(sweep, taper, b , c_r, c_t, c_MAC, dihedral, sweep_LE)
print(c_d0)
print(c_d_0_initial)
print(e_initial)
print(L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D )
print("THIS ONE", y_spanwise, xlemac, lengthMAC)

#____Plots____

plotname = f"Matching_Diagram_iteration_{i}.png"

fig, dx = plt.subplots()
dx.plot(fv.wing_loading_cs, loads_cruise_speed, label="Cruise speed")
dx.plot(fv.wing_loading, loads_climb_grad_119, label="Climb gradient CS25.119")
dx.plot(fv.wing_loading, loads_climb_grad_121a, label="Climb gradient CS25.121a")
dx.plot(fv.wing_loading, loads_climb_grad_121b, label="Climb gradient CS25.121b")
dx.plot(fv.wing_loading, loads_climb_grad_121c, label="Climb gradient CS25.121c")
dx.plot(fv.wing_loading, loads_climb_grad_121d, label="Climb gradient CS25.121d")
dx.plot(fv.wing_loading, loads_to_field, label="Take-off field length")
dx.axvline(loads_minimum_speed[0], color = "gold", label = "Minimum speed")
dx.axvline(loads_landing_field_length[0], color = "black", label = "Landing field length")
dx.set_xlim([0, 9000])
dx.set_ylim([0.1, 0.6])
dx.set_xlabel("W/S - [N/m2]")
dx.set_ylabel("T/W - [N/N]")
dx.set_title("Matching Diagram")
legend = dx.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
plt.savefig(plotname, bbox_extra_artists=(legend,), bbox_inches='tight')

