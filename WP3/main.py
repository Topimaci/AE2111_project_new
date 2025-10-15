import WP2.master_range_mass as mrm
import variables.dynamic_variables as dv
import functions.Class_II_weight_estimations as c2w
import WP2.main as ma
import math as m
import functions.Planform_DESIGN1 as pd
import numpy as np
import functions.Drag_calculations_class_II as D2

mass_to_new = mrm.m_MTO
mass_oe_new = mrm.m_oe
mass_fuel_new = mrm.m_f_des

S_wing = dv.S_w
S_wing_new = S_wing

density_cr = 0.2872
velocity_cr = 200.629


q = c2w.pas_to_psi(0.5 * 0.2872 * fv.v_cr ** 2)
sweep = ma.sweep_true
taper = ma.taper
t_c = ma.thickness_to_chord
AR = ma.AR
b = c2w.m_to_ft()                                               # wing span
S_w = c2w.m2_to_ft2(dv.S_w)                                           # wing area
W_fuel = c2w.kg_to_lb()                                         # Wing fuel weight
N_z = 1.5 * 3.8                                             # Ultimate load factor
N_l =                                                       # Load factor landing
W_des = c2w.kg_to_lb()                                          # Gross design weight
S_wfus = c2w.m2_to_ft2()                                        # Wetted area fuselage
W_l = c2w.kg_to_lb()                                            # Weight landing
M = ma.M_cr                                                 # Mach number
W_uav = c2w.kg_to_lb()                                          # Weight uninstalled avionics
N_pers = 11                                                  # Number personell (Assuming pilot, copilot and one flight attendent)
V_pr = c2w.m3_to_ft3(1.92/2*m.pi*7)                                          # Volume of pressurized section
P_delta = c2w.pas_to_psi()                                      # Cabin pressure
W_press = 11.9 + (V_pr * P_delta) ** 0.271                  # Penalization due to Pressure difference, UNITS?????????????????????????????

while S_wing_new/S_wing <= 0.05:

    

    ##Class II mass calculaitons

    W_wing = c2w.wing_weight(S_w, W_fuel, AR, q, taper, t_c, sweep, N_z, W_des)
    W_htail = c2w.horizontal_tail_weight(N_z, W_des, q, taper, c2w.m2_to_ft2(tail_area), t_c, sweep, htailsweep, h_tailtaper, AR)
    W_vtail = c2w.vertical_tail_weight(N_z, W_des, q, c2w.m2_to_ft2(tail_area), t_c, sweep, sweep_tail, taper_tail, AR)
    W_fuse = c2w.fuselage_weight(S_wfus, N_z, W_des, c2w.m_to_ft(tail_distance), dv.L_over_D_max, q, W_press)
    W_mLG = c2w.main_landing_gear_weight(N_l, W_l, c2w.m_to_in(length_main_gear))
    W_nLG = c2w.nose_landing_gear_weight(N_l, W_l, c2w.m_to_in(length_nose_gear))
    W_eng = c2w.engine_weight(c2w.kg_to_lb(weight_engine), 2)
    W_fs = c2w.fuel_system_weight(c2w.liters_to_gal(tot_fuel_vol), c2w.liters_to_gal(int_tank_vol), number_fueltanks, 2)
    W_fc, W_hyd = c2w.flight_control_and_hydraulics_weight(c2w.m_to_ft(L_fuselage), b, N_z, W_des)
    W_elec, W_avi, W_aircon, W_furn = c2w.electronics_and_avionics_aircondition_furnishings_weight(W_fs, W_uav, W_des, N_pers, M)

    Total_Class_II_Weight = W_wing + W_htail + W_vtail + W_fuse + W_mLG + W_nLG + W_eng + W_fs + W_fc + W_hyd + W_elec + W_avi + W_aircon + W_furn
    ##MTO and OE to be added hellyeah

    ##Planfooooooooooooooooooooooorm calcs


    C_L_des = pd.C_L_design(M_MTO, velocity_cr, density_cr, S_w)
    sweep_LE_DD = pd.sweep_drag_divergence(C_L_des)

    taper, span, chord_root, chord_tip, chord_MAC, dihedral = pd.calculate_geometric_parameters_wing(S_w,AR, M)

    sweep_c4 = pd.sweep_converter(sweep_LE_DD, chord_root, taper, 1/4, span)

    MAC_y, MAC_x = pd.calculate_MAC_position(span, chord_root, chord_tip, sweep_c4)

    ##Drag show

    CD_0_Fus = D2.fuselage_drag_coefficient(S_w, density_cr, velocity_cr, chord_MAC, 1.79*10**(-5), length_fus, diameter_fus, length_cock, length_cyli, length_tail, M, upsweep_tail, Base_are)


