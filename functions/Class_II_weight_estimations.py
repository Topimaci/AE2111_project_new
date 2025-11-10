import math as m
import variables.fixed_values as fv
import variables.dynamic_variables as dv
import functions.Planform_DESIGN1 as Pl
import WP2.main as ma

#--------------Convert Inputs to Imperial Units-------------------------------------------------------------------
def m_to_ft(m):
    ft = 3.28084 * m
    return ft

def m2_to_ft2(m2):
    ft2 = 10.7639 * m2
    return ft2

def kg_to_lb(kg):
    lb = 2.20462 * kg
    return lb

def liters_to_gal(lit):
    gal = 0.264172 * lit
    return gal

def m3_to_ft3(m3):
    ft3 = 35.3147 * m3
    return ft3

def m_to_in(m):
    inch = 39.3701 * m
    return inch

def pas_to_psi(pas):
    psi = 0.000145038 * pas
    return psi

def lb_to_kg(lb):
    kg = 0.45359237 * float(lb)
    return kg



#--------------Class II aircraft weight estimations, formulas from Reymar-----------------------------------------
def wing_weight(wing_area, wing_fuel_weight, aspect_ratio, dynamic_pressure, taper_ratio, tc_ratio, quarter_sweep, ultimate_load, design_gross_weight):
    W_wing = 0.036 * wing_area ** 0.758 * wing_fuel_weight ** 0.0035 * (aspect_ratio / (m.cos(quarter_sweep)) ** 2) ** 0.6 * dynamic_pressure ** 0.006 * taper_ratio ** 0.04 * (100 * tc_ratio / m.cos(quarter_sweep)) ** (-0.3) * (ultimate_load * design_gross_weight) ** 0.49
    return W_wing


def horizontal_tail_weight(ultimate_load, design_gross_weight, dynamic_pressure, h_tail_area, tc_ratio, quarter_sweep, htail_quarter_sweep, h_tail_taper, aspect_ratio):
    htail_weight = 0.016 * (ultimate_load * design_gross_weight) ** 0.414 * dynamic_pressure ** 0.168 * h_tail_area ** 0.896 * (100 *  tc_ratio / (m.cos(quarter_sweep)) ** 2) ** (-0.12) * (aspect_ratio / (m.cos(htail_quarter_sweep)) ** 2) ** 0.043 * h_tail_taper ** (-0.02)
    return htail_weight

def vertical_tail_weight(ultimate_load, design_gross_weight, dynamic_pressure, v_tail_area, tc_ratio, quarter_sweep, vtail_quarter_sweep, v_tail_taper, aspect_ratio, H_t_H_v_ratio):
    vtail_weight = 0.073 * (1 + 0.2 * H_t_H_v_ratio) * (ultimate_load * design_gross_weight) ** 0.376 * dynamic_pressure ** 0.122 * v_tail_area ** 0.873 * (100 * tc_ratio / m.cos(quarter_sweep)) ** (-0.49) * (aspect_ratio / (m.cos(vtail_quarter_sweep)) ** 2) ** 0.357 * v_tail_taper ** 0.039
    return vtail_weight

def fuselage_weight(fuselage_wetted_A, ultimate_load, design_gross_weight, tail_distance, L_fuselage, D_fuselage, dynamic_pressure, W_pressure):        # D_fuselage is the structural depth of the fuselage
    fuse_weight = 0.052 * (fuselage_wetted_A) ** 1.086 * (ultimate_load * design_gross_weight) ** 0.177 * tail_distance ** (-0.051) * L_fuselage / D_fuselage ** (-0.072) * dynamic_pressure ** 0.241 + W_pressure
    return fuse_weight

def main_landing_gear_weight(ultimate_landing_load, weight_landing, length_main_gear):
    main_LG_weight = 0.095 * (ultimate_landing_load * weight_landing) ** 0.768 * (length_main_gear / 12) ** 0.409
    return main_LG_weight

def nose_landing_gear_weight(ultimate_landing_load, weight_landing, length_nose_gear):
    nose_LG_weight = 0.125 * (ultimate_landing_load * weight_landing) ** 0.566 * (length_nose_gear / 12) ** 0.845
    return nose_LG_weight

def engine_weight(weight_1engine, number_engines):
    engine_weight = 2.575 * weight_1engine ** 0.922 * number_engines
    return engine_weight

def fuel_system_weight(tot_fuel_vol, int_tank_vol, number_fueltanks, number_engines):
    fuel_sys_weight = 2.49 * tot_fuel_vol ** 0.726 * (1 / (1 + int_tank_vol / tot_fuel_vol)) ** 0.363 * number_fueltanks ** 0.242 * number_engines ** 0.157
    return fuel_sys_weight

def flight_control_and_hydraulics_weight(L_fuselage, span_b, ultimate_load, design_gross_weight):
    fligth_controls_weight = 0.053 * L_fuselage ** 1.536 * span_b ** 0.371 * (ultimate_load * design_gross_weight * 10**(-4)) ** 0.8
    hydraulics_weight = 0.001 * design_gross_weight
    return fligth_controls_weight, hydraulics_weight

def electronics_and_avionics_aircondition_furnishings_weight(weight_fuel_system, uninstalled_avionics_W, design_gross_weight, number_personell, Mach):
    avionics_weight = 2.117 * uninstalled_avionics_W ** 0.933
    electronics_weight = 12.57 * (avionics_weight + weight_fuel_system) ** 0.51 
    aircondition_weight = 0.265 * design_gross_weight ** 0.52 * number_personell ** 0.68 * avionics_weight ** 0.17 * Mach ** 0.08
    furnishings_weight = 0.0582 * design_gross_weight - 65
    return electronics_weight, avionics_weight, aircondition_weight, furnishings_weight






def max_takeoff_mass(OEM, crew, W_fuel, W_payload):         # Crew is equal to number of crew 
    MTO = OEM + crew * 90 * 1.05 + W_fuel + W_payload       # in kg!!!!!!!!!!!!!!!!!!
    return MTO
