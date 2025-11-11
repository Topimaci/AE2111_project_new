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
import functions.Range_calculations as range

#from scipy.optimize import fsolve

MTOW = 12520
OEW = 7608
W_fuel = 4161
mass_landing = MTOW-W_fuel
S_wing = 49.63

Running = True
#change dv.S_w to S_wing in all code DONE

print("Fuel:" , W_fuel)
density_cr = 0.2872
velocity_cr = 200.629

visc = 1.79*(10**(-5))

##FUselage variables
length_fus = 15.2
diameter_fus = 2.2
length_cock = 3
length_cyli = 10
length_tail = 2.2

#
number_fueltanks = 3



q = c2w.pas_to_psi(0.5 * 0.2872 * fv.v_cr ** 2)
sweep = 11.852

taper = ma.taper
t_c = ma.thickness_to_chord
AR = ma.AR
span = ma.b                                              # wing span
N_z = 1.5 * 3.8                                             # Ultimate load factor
N_l = 2.5                                                      # Load factor landing
W_des = MTOW-0.5*W_fuel                                          # Gross design weight
S_wfus = (m.pi*diameter_fus/4)*(1/(3*length_cock**2)*((4*length_cock**2 + (diameter_fus**2)/4)**1.5 -(diameter_fus**3 /8))-diameter_fus+4*length_cyli+2*m.sqrt(length_tail**2 + (diameter_fus**2)/4))                                        # Wetted area fuselage

M = ma.M_cr                                                 # Mach number
N_pers = 11                                                  # Number personell (Assuming pilot, copilot and one flight attendent)
V_pr = c2w.m3_to_ft3(1.92/2*m.pi*7)                                          # Volume of pressurized section
P_delta = 8                                      # Cabin pressure
W_press = 11.9 + (V_pr * P_delta) ** 0.271                  # Penalization due to Pressure difference, UNITS?????????????????????????????
chord_MAC = 2.39
chord_root = 3.269

fuel_mass_fraction = 0.333
m_OE = 0.607
m_wing = 0.135
m_fus = 0.105
m_t = 0.043
m_eng = 0.1
m_nac = 0.018
m_lg = 0.036
m_fe = 0.17
m_unacc = 0

AR_h = 3.5                                                      #Aspect ratio horizontal tail
AR_v = 1.5                                                      #Aspect ratio vertical tail
h_tailtaper = 0.7                                               #Horizontal tail taper ratio
v_tailtaper = 0.9                                               #Horizontal tail taper ratio 
v_tailsweep = 20                                                #Vertical Tail sweep, degrees
h_tailsweep = 15                                                #Horizontal Tail sweep, degrees

t_w = dv.designtw
w_s = dv.designws


sweep_t_c_max = pd.sweep_converter(sweep, chord_root, taper, 0.3, span)
i = 0

while Running == True:
    i = i+1

    tail_area_v, tail_area_h, tail_distance =  es.calculate_tail_surface_areas(S_wing, span, chord_MAC, chord_root, MTOW, fuel_mass_fraction, m_OE, m_wing, m_fus, m_t, m_eng, m_nac, m_fe)
    



    ##Class II mass calculaitons

    #EEEEEEEEEEEEEEEEEEEEEEEEEngine

    engine = et.engine_required(t_w*9.81*MTOW)
    print("engine", engine)

    sweep_c4 = pd.sweep_converter(sweep, chord_root, taper, 0.25, span)


    W_wing = c2w.lb_to_kg(c2w.wing_weight(c2w.m2_to_ft2(S_wing),  c2w.kg_to_lb(W_fuel), AR, c2w.pas_to_psi(q), taper, t_c, sweep, N_z,  c2w.kg_to_lb(W_des)))
    W_htail = c2w.lb_to_kg(c2w.horizontal_tail_weight(N_z,  c2w.kg_to_lb(W_des), c2w.pas_to_psi(q), c2w.m2_to_ft2(tail_area_h), t_c, sweep, h_tailsweep, h_tailtaper, AR_h))
    W_vtail = c2w.lb_to_kg(c2w.vertical_tail_weight(N_z,  c2w.kg_to_lb(W_des), c2w.pas_to_psi(q), c2w.m2_to_ft2(tail_area_v), t_c, sweep, v_tailsweep, v_tailtaper, AR_v, 0))
    W_fuse = c2w.lb_to_kg(c2w.fuselage_weight(c2w.m2_to_ft2(S_wfus), N_z,  c2w.kg_to_lb(W_des), c2w.m_to_ft(tail_distance), c2w.m_to_ft(length_fus), c2w.m_to_ft(diameter_fus), c2w.pas_to_psi(q),  c2w.kg_to_lb(W_press)))
    W_mLG = c2w.lb_to_kg(c2w.main_landing_gear_weight(N_l,  c2w.kg_to_lb(mass_landing), c2w.m_to_in(0.74)))  ###assumed 1.5
    W_nLG = c2w.lb_to_kg(c2w.nose_landing_gear_weight(N_l,  c2w.kg_to_lb(mass_landing), c2w.m_to_in(0.74))) ###assumed 1.2
    W_eng = c2w.lb_to_kg(c2w.engine_weight(c2w.kg_to_lb(engine[2]), 2))
    W_fs = c2w.lb_to_kg(c2w.fuel_system_weight(c2w.liters_to_gal(W_fuel/800*1000), c2w.liters_to_gal(W_fuel/800*1000), number_fueltanks, 2))
    W_fc, W_hyd = c2w.flight_control_and_hydraulics_weight(c2w.m_to_ft(length_fus), c2w.m_to_ft(span), N_z,  c2w.kg_to_lb(W_des))
    W_fc = c2w.lb_to_kg(W_fc)
    W_hyd = c2w.lb_to_kg(W_hyd)
    W_elec, W_avi, W_aircon, W_furn = c2w.electronics_and_avionics_aircondition_furnishings_weight( c2w.kg_to_lb(W_fs),  c2w.kg_to_lb(607),  c2w.kg_to_lb(W_des), N_pers, M)
    W_elec = c2w.lb_to_kg(W_elec)
    W_avi = c2w.lb_to_kg(W_avi)
    W_aircon = c2w.lb_to_kg(W_aircon)
    W_furn = c2w.lb_to_kg(W_furn)
    W_payload  = 1010

    OEW = W_wing + W_htail + W_vtail + W_fuse + W_mLG + W_nLG + W_eng + W_fs + W_fc + W_hyd + W_elec + W_avi + W_aircon + W_furn

    MTOW = c2w.max_takeoff_mass(OEW, 3, W_fuel, W_payload)



    ##MTO and OE to be added hellyeah

    ##Planfooooooooooooooooooooooorm calcs


    C_L_des = pd.C_L_design(MTOW, W_fuel, velocity_cr, density_cr, S_wing)
    #sweep_LE_DD = pd.sweep_drag_divergence(C_L_des)  
    #print("sweep:" ,sweep_LE_DD)
    ##sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE

    sweep_false, taper, span, chord_root, chord_tip, chord_MAC, dihedral, sweep_LE_false = pd.calculate_geometric_parameters_wing(S_wing,AR, M)
    sweep_htail_false, htail_taper, htail_span, htail_chord_root, htail_chord_tip, htail_chord_MAC, htail_dihedral, sweep_LE_h_false = pd.calculate_geometric_parameters_wing(tail_area_h,3.5,M)
    sweep_vtail_false, vtail_taper, vtail_span, vtail_chord_root, vtail_chord_tip, vtail_chord_MAC, vtail_dihedral, sweep_LE_v_false = pd.calculate_geometric_parameters_wing(tail_area_v,1.5,M)



    cf, S_flap = hld.HLD(S_wing, sweep, span, chord_tip, chord_root)



    print(f"[HLD] cf/c={cf:.3f}, S_flap={S_flap:.3f} m² (S_wing={S_wing:.1f} m²)")

    #y_spanwise, xlemac, lengthMAC
    ##Drag show

    ####assumptions for base area and upsweep
    upsweep_tail = 0.09 ##20deg
    Base_area = 0.0314 

    CD_0_Fus = D2.fuselage_drag_coefficient(S_wing, density_cr, velocity_cr, chord_MAC, visc, length_fus, diameter_fus, length_cock, length_cyli, length_tail, M, upsweep_tail, Base_area)
    CD_0_Wing = D2.wing_drag_coefficient(0.14,0.30,pd.sweep_converter(sweep_t_c_max,chord_root, taper, 0.30, span),S_wing, density_cr, velocity_cr, chord_MAC, visc, M)
    CD_0_Htail = D2.horizontal_tail_drag_coefficient(S_wing, 0.12,0.3, pd.sweep_converter(25, htail_chord_root, htail_taper, 0.3, htail_span), tail_area_h, density_cr, velocity_cr, htail_chord_MAC, visc, M)
    CD_0_Vtail = D2.vertical_tail_drag_coefficient(S_wing, 0.12, 0.3, pd.sweep_converter(20, vtail_chord_root, vtail_taper, 0.3, vtail_span),tail_area_v, density_cr, velocity_cr, vtail_chord_MAC, visc, M)
    CD_0_Nacelle = D2.nacelle_drag_coefficient(S_wing, density_cr, velocity_cr, visc, engine[5], engine[4], M)
    CD_0_surf = CD_0_Fus + CD_0_Wing + 1.06*(CD_0_Htail + CD_0_Vtail) + CD_0_Nacelle
                                ### 1.06 = IFc
    ## Have to decide which ones count to fwhich configuration, also CD_wheelwell times 3 or only 1 or what??
    #CD_Wheelwell = D2.C_D_landing_gear_whells(fuselage_height, width_tire_and_strut, height_strut, width_strut, height_gear, width_gear, S_wing)
    CD_flap = D2.flap_drag_coefficient(cf/chord_MAC,S_flap,S_wing, 40)

    CD_0_misc = CD_flap
    sweep_half = pd.sweep_converter(sweep, chord_root, taper, 0.5, span)
    CD_ind_clean, e_clean, AR_new_clean = D2.induced_drag(AR, sweep_half, C_L_des, 0, 1.3, span, 10000000000000) ########find reference wignlet
    CD_ind_Landing, e, AR_new = D2.induced_drag(AR, sweep_half, 2.59, 40, 1.3, span, 100000000000000000000000)
                                                        ### 2.59 assumed from WP2
    CD_wave = D2.wave_C_D(M, 0.70)

    CD_0_final = CD_0_surf+CD_0_misc+0.03*(CD_0_surf+CD_0_misc)

    CD_total_clean = CD_0_surf*1.03+CD_ind_clean+CD_wave
    CD_total_landing = CD_0_final + CD_ind_Landing
    print("Cdtotalclen" , CD_total_clean)
    print(f"""
    ==================== DRAG COMPONENT SUMMARY ====================

    Parasite Drag Components:
    CD₀ (Fuselage):        {CD_0_Fus:.6f}
    CD₀ (Wing):            {CD_0_Wing:.6f}
    CD₀ (H. Tail):         {CD_0_Htail:.6f}
    CD₀ (V. Tail):         {CD_0_Vtail:.6f}
    CD₀ (Nacelle):         {CD_0_Nacelle:.6f}
    ---------------------------------------------------------------
    CD₀ (Surface total):   {CD_0_surf:.6f}
    CD₀ (Flaps):           {CD_flap:.6f}
    CD₀ (Misc.):           {CD_0_misc:.6f}
    CD₀ (Final):           {CD_0_final:.6f}

    Induced and Wave Drag:
    CDᵢ (Clean):           {CD_ind_clean:.6f}
    CDᵢ (Landing):         {CD_ind_Landing:.6f}
    Oswald e (Clean):      {e_clean:.4f}
    Oswald e (Landing):    {e:.4f}
    AR (Clean):            {AR_new_clean:.3f}
    AR (Landing):          {AR_new:.3f}
    CD_wave:               {CD_wave:.6f}

    Total Drag:
    CD_total (Clean):      {CD_total_clean:.6f}
    CD_total (Landing):    {CD_total_landing:.6f}

    ===============================================================
    """)

    #LIFT OVER DRAG
    L_over_D = C_L_des/CD_total_clean
    print("LoverD", L_over_D)

    #Fuel mass fraction
    print("CL design", C_L_des)
 
    R_lost = range.R_lost_function(L_over_D, fv.h_cr, velocity_cr)
    R_eq_res = range.R_eq_res_function(fv.R_div, fv.t_E, velocity_cr)
    R_eq = range.R_eq_function(fv.R_des, R_lost, fv.f_cont, R_eq_res)
    eta_j = range.eta_j_function(fv.e_f, velocity_cr, engine[8])
    fuel_mass_fraction = range.fuel_mass_fraction_function(R_eq, eta_j, fv.e_f, L_over_D)
    print("fmf", fuel_mass_fraction)
    W_fuel = fuel_mass_fraction * MTOW

    mass_landing = MTOW - W_fuel

    m_OE = OEW / MTOW
    m_wing = W_wing / MTOW
    m_fus = W_fuse / MTOW
    m_t = (W_vtail + W_htail)/MTOW
    m_eng = (W_eng)/ MTOW
    m_nac = (engine[3]-engine[2])/MTOW
    m_fe = 0.17
    m_unacc = 1-(m_OE+m_wing+m_fus+m_t+m_eng+m_nac+m_fe)

    Total_volume_wing, percentage_fuel_in_wing, fuel_mass_in_wing, Total_volume_fuel_needed = fuelv.fuel_volume(chord_root, taper, span, W_fuel)
    print(f"""
    Total Fuel Volume: {Total_volume_fuel_needed}
    Total Fuel Volume in wing: {Total_volume_wing}
    Percentage of Fuel in wing: {percentage_fuel_in_wing}
    Fuel mass in wing: {fuel_mass_in_wing}        
    """)
    #matching diagram


    loads_minimum_speed = ms.Minimum_speed(1.225, engine[7], 66, 2.59) ### 2.59 is CL_max
    loads_landing_field_length = lfl.landing_field_length(mass_landing/MTOW, 700, 1.225, 2.59)  ### changed landing mass fraction from constant to iterative
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

    intersection_tola = np.interp(loads_landing_field_length, fv.wing_loading, loads_to_field)

    from scipy.optimize import fsolve

    def f(w):
        return np.interp(w, fv.wing_loading, loads_to_field) - np.interp(w, fv.wing_loading, loads_cruise_speed)

    try:
        intersection_tocr = fsolve(f, x0=fv.wing_loading.mean())[0]
    except Exception:
        intersection_tocr = np.nan

    wingload_1 = intersection_tocr

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

    print(f"""
    Span:        {span:.3f} m
    Chord root:  {chord_root:.3f} m
    Taper:       {taper:.3f}
    LE Sweep:    {sweep:.2f}°
    Wing Area:   {S_wing_new:.3f} m²
    CL Design:   {C_L_des:.3f}
    MTOW:        {MTOW:.1f} kg
    Fuel mass:   {W_fuel:.1f} kg
    """)


    if abs(S_wing_new - S_wing)/S_wing < 0.0000000000001 or i == 100:

        print(i)
        Running = False

    else: S_wing = S_wing_new