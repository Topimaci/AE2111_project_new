import WP2.master_range_mass as mrm
import variables.dynamic_variables as dv
import functions.Class_II_weight_estimations as c2w
import WP2.main as ma
import math as m
import functions.Planform_DESIGN1 as pd
import numpy as np
import functions.Drag_calculations_class_II as D2
import functions.Empennage_sizing as es
import functions.Engine_types as et
import WP2.HLDs as hld
import functions.Minimum_speed as ms
import functions.Climb_rate as cr
import functions.Cruise_speed as cs
import functions.Landing_field_length as lfl
import functions.climb_grad as cg
import functions.Take_off_distance as td
import variables.fixed_values as fv
import functions.Fuel_Volume as fuelv

mass_to_new = mrm.m_MTO
mass_oe_new = mrm.m_oe
mass_fuel_new = mrm.m_f_des

S_wing = dv.S_w

Running = True
#change dv.S_w to S_wing in all code DONE


density_cr = 0.2872
velocity_cr = 200.629

visc = 1.79*10**(-5)

##FUselage variables
length_fus = 15.2
diameter_fus = 2.2
length_cock = 3
length_cyli = 10
length_tail = 2.2

#
number_fueltanks = 3



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

AR_h = 3.5                                                      #Aspect ratio horizontal tail
AR_v = 1.5                                                      #Aspect ratio vertical tail
h_tailtaper = 0.7                                               #Horizontal tail taper ratio
v_tailtaper = 0.9                                               #Horizontal tail taper ratio 
v_tailsweep = 20                                                #Vertical Tail sweep, degrees
h_tailsweep = 15                                                #Horizontal Tail sweep, degrees

t_w = dv.designtw
w_s = dv.designws

while Running == True:

    tail_area_v, tail_area_h, tail_distance =  es.calculate_tail_surface_areas(S_w, span, chord_MAC, chord_root, MTOW, fuel_mass_fraction, m_OE, m_wing, m_fus, m_t, m_eng, m_nac, m_lg, m_fe, m_unacc)
    



    ##Class II mass calculaitons

    #EEEEEEEEEEEEEEEEEEEEEEEEEngine

    engine = et.engine_required(t_w*9.81*mass_to_new)

    Total_volume_wing, percentage_fuel_in_wing, fuel_mass_in_wing, Total_volume_fuel_needed = fuelv.fuel_volume


    W_wing = c2w.wing_weight(S_w, W_fuel, AR, q, taper, t_c, sweep, N_z, W_des)
    W_htail = c2w.horizontal_tail_weight(N_z, W_des, q, taper, c2w.m2_to_ft2(tail_area_h), t_c, sweep, h_tailsweep, h_tailtaper, AR_h)
    W_vtail = c2w.vertical_tail_weight(N_z, W_des, q, c2w.m2_to_ft2(tail_area_v), t_c, sweep, v_tailsweep, v_tailtaper, AR_v)
    W_fuse = c2w.fuselage_weight(S_wfus, N_z, W_des, c2w.m_to_ft(tail_distance), dv.L_over_D_max, q, W_press)
    W_mLG = c2w.main_landing_gear_weight(N_l, W_l, c2w.m_to_in(length_main_gear))
    W_nLG = c2w.nose_landing_gear_weight(N_l, W_l, c2w.m_to_in(length_nose_gear))
    W_eng = c2w.engine_weight(c2w.kg_to_lb(engine[7]), 2)
    W_fs = c2w.fuel_system_weight(c2w.liters_to_gal(Total_volume_fuel_needed*1000), c2w.liters_to_gal(Total_volume_fuel_needed*1000), number_fueltanks, 2)
    W_fc, W_hyd = c2w.flight_control_and_hydraulics_weight(c2w.m_to_ft(length_fus), b, N_z, W_des)
    W_elec, W_avi, W_aircon, W_furn = c2w.electronics_and_avionics_aircondition_furnishings_weight(W_fs, W_uav, W_des, N_pers, M)
    W_payload  = 1010

    OEW = W_wing + W_htail + W_vtail + W_fuse + W_mLG + W_nLG + W_eng + W_fs + W_fc + W_hyd + W_elec + W_avi + W_aircon + W_furn

    MTOW = c2w.max_takeoff_mass(OEW, 3, W_fuel, W_payload)



    fuel_mass_fraction = 
    m_OE = OEW / MTOW
    m_wing = W_wing / MTOW
    m_fus = W_fuse / MTOW
    m_t = (W_vtail + W_htail)/MTOW
    m_eng = (W_eng)/ MTOW
    m_nac = 
    m_lg = 
    m_fe = 
    m_unacc = 

    ##MTO and OE to be added hellyeah

    ##Planfooooooooooooooooooooooorm calcs


    C_L_des = pd.C_L_design(MTOW, mass_fuel_new, velocity_cr, density_cr, S_wing)
    C_L_land = pd.C_L_design(M)
    sweep_LE_DD = pd.sweep_drag_divergence(C_L_des)

    taper, span, chord_root, chord_tip, chord_MAC, dihedral = pd.calculate_geometric_parameters_wing(S_wing,AR, M)
    htail_taper, htail_span, htail_chord_root, htail_chord_tip, htail_chord_MAC, htail_dihedral = pd.calculate_geometric_parameters_wing(es.S_h,3.5,M)
    vtail_taper, vtail_span, vtail_chord_root, vtail_chord_tip, vtail_chord_MAC, vtail_dihedral = pd.calculate_geometric_parameters_wing(es.S_v,3.5,M)

    sweep_c4 = pd.sweep_converter(sweep_LE_DD, chord_root, taper, 1/4, span)

    MAC_y, MAC_x = pd.calculate_MAC_position(span, chord_root, chord_tip, sweep_c4)
    cf, S_flap = hld.HLD(S_wing_new, sweep_LE_DD, span, chord_tip, chord_root)

    ##Drag show

    CD_0_Fus = D2.fuselage_drag_coefficient(S_wing, density_cr, velocity_cr, chord_MAC, visc, length_fus, diameter_fus, length_cock, length_cyli, length_tail, M, upsweep_tail, Base_area)
    CD_0_Wing = D2.wing_drag_coefficient(0.14,0.378,pd.sweep_converter(sweep_LE_DD,chord_root, taper, 0.378, span),S_wing, density_cr, velocity_cr, chord_MAC, visc, M)
    CD_0_Htail = D2.horizontal_tail_drag_coefficient(0.12,0.3, pd.sweep_converter(25, htail_chord_root, htail_taper, 0.3, htail_span), es.S_h, density_cr, velocity_cr, htail_chord_MAC, visc, M)
    CD_0_Vtail = D2.vertical_tail_drag_coefficient(0.12, 0.3, pd.sweep_converter(20, vtail_chord_root, vtail_taper, 0.3, vtail_span),es.S_v, density_cr, velocity_cr, vtail_chord_MAC, visc, M)
    CD_0_Nacelle = D2.nacelle_drag_coefficient(S_wing, density_cr, velocity_cr, visc, engine[5], engine[4], M)
    CD_0_surf = CD_0_Fus + CD_0_Wing + CD_0_Htail + CD_0_Vtail + CD_0_Nacelle

    ## Have to decide which ones count to fwhich configuration, also CD_wheelwell times 3 or only 1 or what??
    CD_Wheelwell = D2.C_D_landing_gear_whells(fuselage_height, width_tire_and_strut, height_strut, width_strut, height_gear, width_gear, S_wing)
    CD_flap = D2.flap_drag_coefficient(cf/chord_MAC,S_flap,S_wing_new, 40)

    CD_0_misc = CD_Wheelwell + CD_flap

    CD_ind_clean, e, AR_new = D2.induced_drag(AR, sweep_LE_DD, C_L_des, 0, 1, span, 10000000000000)
    CD_ind_Landing, e, AR_new = D2.induced_drag(AR, sweep_LE_DD, C_L_land, 40, 1, span, 100000000000000000000000)

    CD_wave = D2.wave_C_D(M, 0.68)

    CD_0_final = CD_0_surf+CD_0_misc+0.03*(CD_0_surf+CD_0_misc)

    CD_total_clean = CD_0_surf*1.03+CD_ind_clean+CD_wave
    CD_total_landing = CD_0_final + CD_ind_Landing

    #matching diagram

    loads_minimum_speed = ms.Minimum_speed(1.225, engine[7], 66, C_L_land)
    loads_landing_field_length = lfl.landing_field_length(mass_landing/mass_to_new, 700, 1.225, C_L_land)
    loads_cruise_speed = cs.cruise_speed(0.95,0.24,fv.wing_loading_cs,CD_0_surf, 0.2872, velocity_cr, AR_new, e)
    loads_climb_rate = cr.climb_rate(fv.wing_loading)
    loads_climb_grad_119 = cg.climb_grad(fv.wing_loading, fv.mass_fraction_119, fv.cg_119, fv.C_d0_119, fv.e_119, fv.AR, 1.225, fv.C_l_at_max_climb_gradient_119, fv.B)
    loads_climb_grad_121a = 2*cg.climb_grad(fv.wing_loading, fv.mass_fraction_121a, fv.cg_121a,fv.C_d0_121a, fv.e_121a, fv.AR, 1.225, fv.C_l_at_max_climb_gradient_121a, fv.B)
    loads_climb_grad_121b = 2*cg.climb_grad(fv.wing_loading, fv.mass_fraction_121b, fv.cg_121b, fv.C_d0_121b, fv.e_121b, fv.AR, 1.225, fv.C_l_at_max_climb_gradient_121b, fv.B)
    loads_climb_grad_121c = 2*cg.climb_grad(fv.wing_loading, fv.mass_fraction_121c, fv.cg_121c, fv.C_d0_121c, fv.e_121c, fv.AR, 1.225, fv.C_l_at_max_climb_gradient_121c, fv.B)
    loads_climb_grad_121d = 2*cg.climb_grad(fv.wing_loading, fv.mass_fraction_121d, fv.cg_121d, fv.C_d0_121d, fv.e_121d, fv.AR, 1.225, fv.C_l_at_max_climb_gradient_121d, fv.B)
    V2 = td.find_V_2(fv.wing_loading, fv.density_takeoff, fv.CL_takeoff)
    Mach = td.find_Mach(V2, fv.temp_takeoff)
    theta, delta = td.find_theta_delta(fv.temp_takeoff, Mach)
    alpha = td.find_alpha_t(delta, Mach, fv.B)
    loads_to_field = td.take_off_distance(alpha, fv.wing_loading, fv.takeoff_field, fv.density_takeoff, fv.oswald_efficiency, fv.AR)

    intersection_tocr = np.intersect1d(loads_to_field, loads_cruise_speed)
    intersection_tola = loads_to_field[np.where(fv.wing_loading == loads_landing_field_length)]

    wingload_1 = fv.wing_loading[np.where(loads_to_field == intersection_tocr)]
    wingload_2 = loads_landing_field_length

    if wingload_1*1.05 <= wingload_2 and wingload_1 <= wingload_2:

        w_s_new = wingload_2*0.95
        t_w_new = intersection_tola * 1.05
    
    elif wingload_1*1.05 >= wingload_2 and wingload_1 <= wingload_2: 

        w_s_new = wingload_1
        t_w_new = intersection_tocr * 1.05
    else:

        if intersection_tocr >= intersection_tola:

            w_s_new = wingload_2*0.95
            t_w_new = intersection_tocr*1.05

        else:

            w_s_new = wingload_2*0.95
            t_w_new = intersection_tola*1.05

    W_to = MTOW * 9.81
    S_wing_new = W_to/w_s_new

    if S_wing_new/S_wing <= 0.05 or S_wing_new/S_wing >= 0.05:

        Running = False

    else: S_wing = S_wing_new