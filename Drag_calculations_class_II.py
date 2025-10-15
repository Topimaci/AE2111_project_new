import math
from math import pi as pi

#25% of fuselage is laminar rest is turbulent, from table
def fuselage_drag_coefficient(wing_area, density, velocity, MAC, dynamic_viscosity, l_fuselage, diameter_fuselage, length_cockpit, length_cylinder_part, length_tail_part, Mach, upsweep, Base_area):
    Re_1 = density*velocity*MAC/dynamic_viscosity
    Re_check = 38.21*(l_fuselage/(0.052*10**(-5)))**1.053 #assumed smooth molded composite
    if Re <= Re_check:
        Re = Re_1
    else:
        Re = Re_check
    #laminar friciton coefficient calculations
    C_f_laminar = 1.328/math.sqrt(Re)
    C_f_turbulent = 0.455/((math.log10(Re))**(2.58) * (1+0.144*Mach**2)**(0.65))

    c_f_total = 0.25*C_f_laminar + 0.75*C_f_turbulent

    f = l_fuselage/diameter_fuselage
    FF = (1+60/f**3 + f/400)
    S_wet_fus = (math.pi*diameter_fuselage/4)*(1/(3*length_cockpit**2)*((4*length_cockpit**2 + (diameter_fuselage**2)/4)**1.5 -(diameter_fuselage**3 /8))-diameter_fuselage+4*length_cylinder_part+2*math.sqrt(length_tail_part**2 + (diameter_fuselage**2)/4))
    C_D_0_fuselage_friction = FF*c_f_total*S_wet_fus/wing_area

     #### miscellaneous drag fuselage
    #upsweep                                                          
    C_D_upsweep = (3.83*upsweep**(2.5)*((math.pi*(diameter_fuselage/2)**2)))/wing_area
    C_D_base_drag = ((0.139+0.419*(Mach-0.161)**2)*Base_area)
    C_D_0_fuselage = C_D_0_fuselage_friction + C_D_base_drag + C_D_upsweep
    return C_D_0_fuselage

##35% of wing is laminar
def wing_drag_coefficient(t_over_c, x_over_c_max, sweep_max_t_c, wing_area, density, velocity, MAC, dynamic_viscosity, Mach):
    Re_1 = density*velocity*MAC/dynamic_viscosity
    Re_check = 38.21*(MAC/(0.152*10**(-5)))**1.053 #assumed polished sheet metal
    if Re <= Re_check:
        Re = Re_1
    else:
        Re = Re_check

    C_f_laminar = 1.328/math.sqrt(Re)
    C_f_turbulent = 0.455/((math.log10(Re))**(2.58) * (1+0.144*Mach**2)**(0.65))

    c_f_total = 0.35*C_f_laminar+0.65*C_f_turbulent

    FF = ((1+0.6*t_over_c/x_over_c_max +100*t_over_c**4)*(1.34*Mach**0.18 *(math.cos(sweep_max_t_c))**0.28))
    S_wet_wing = 1.07*2* wing_area 
    IF_c = 1.25 ##interference when connecting wing to fuselage
    C_D_0_wing = FF * c_f_total* IF_c * S_wet_wing

    return C_D_0_wing

#35% is laminar
def horizontal_tail_drag_coefficient(t_over_c_htail, x_over_c_max_htail, sweep_max_t_c_htail, horizontal_tail_area, density, velocity, MAC_htail, dynamic_viscosity, Mach):
    Re_1 = density*velocity*MAC_htail/dynamic_viscosity
    Re_check = 38.21*(MAC_htail/(0.152*10**(-5)))**1.053 #assumed polished sheet metal
    if Re <= Re_check:
        Re = Re_1
    else:
        Re = Re_check

    C_f_laminar = 1.328/math.sqrt(Re)
    C_f_turbulent = 0.455/((math.log10(Re))**(2.58) * (1+0.144*Mach**2)**(0.65))

    c_f_total = 0.35*C_f_laminar+0.65*C_f_turbulent

    FF = ((1+0.6*t_over_c_htail/x_over_c_max_htail +100*t_over_c_htail**4)*(1.34*Mach**0.18 *(math.cos(sweep_max_t_c_htail))**0.28))
    S_wet_wing = 1.07*2* horizontal_tail_area
    IF_c = 1.044 ##interference when connecting wing to fuselage
    C_D_0_horizontal_tail = FF * IF_c* c_f_total * S_wet_wing
    return C_D_0_horizontal_tail


def vertical_tail_drag_coefficient(t_over_c_vtail, x_over_c_max_vtail, sweep_max_t_c_vtail, vertical_tail_area, density, velocity, MAC_vtail, dynamic_viscosity, Mach):
    Re_1 = density*velocity*MAC_vtail/dynamic_viscosity
    Re_check = 38.21*(MAC_vtail/(0.152*10**(-5)))**1.053 #assumed polished sheet metal
    if Re <= Re_check:
        Re = Re_1
    else:
        Re = Re_check

    C_f_laminar = 1.328/math.sqrt(Re)
    C_f_turbulent = 0.455/((math.log10(Re))**(2.58) * (1+0.144*Mach**2)**(0.65))

    c_f_total = 0.35*C_f_laminar+0.65*C_f_turbulent

    FF = ((1+0.6*t_over_c_vtail/x_over_c_max_vtail +100*t_over_c_vtail**4)*(1.34*Mach**0.18 *(math.cos(sweep_max_t_c_vtail))**0.28))
    S_wet_wing = 1.07*2* vertical_tail_area
    IF_c = 1.044 ##interference when connecting wing to fuselage
    C_D_0_vertical_tail = FF * IF_c * c_f_total * S_wet_wing
    return C_D_0_vertical_tail


def nacelle_drag_coefficient(wing_area, density, velocity, dynamic_viscosity, length_nacelle, diameter_nacelle,Mach):
    Re_1 = density*velocity*length_nacelle/dynamic_viscosity
    Re_check = 38.21*(length_nacelle/(0.052*10**(-5)))**1.053 #assumed smooth molded composite
    if Re <= Re_check:
        Re = Re_1
    else:
        Re = Re_check
    #laminar friciton coefficient calculations
    C_f_laminar = 1.328/math.sqrt(Re)
    C_f_turbulent = 0.455/((math.log10(Re))**(2.58) * (1+0.144*Mach**2)**(0.65))

    c_f_total = 0.1*C_f_laminar + 0.9*C_f_turbulent
    ####assumption that only 10% is laminar of nacelle !!!!!!!!!! #####################


    f = length_nacelle/diameter_nacelle
    FF = 1 + 0.35/f
    IF_c = 1.5 
    S_wet_nacelle = math.pi*(diameter_nacelle/2)**2 *length_nacelle
    C_D_0_nacelle = FF*c_f_total* IF_c*S_wet_nacelle/wing_area
    return C_D_0_nacelle

def C_D_landing_gear_whells(h_from_fuselage, width_tire_and_strut, height_strut, width_strut, height_gear, width_gear, wing_area):
    ##assume closed wheels wells
    C_D_s = 0.04955*math.exp(5.615*(width_strut*height_strut+width_gear*height_gear)/(width_tire_and_strut*h_from_fuselage))
    C_D_landing_gear = C_D_s * (h_from_fuselage*width_tire_and_strut)/wing_area
    return C_D_landing_gear
    
def flap_drag_coefficient(c_f_over_c, S_flap, wing_area, deflection_flap):
    C_D_Flap = 0.0074 * c_f_over_c * S_flap/wing_area *(deflection_flap-10)
    return C_D_Flap


# Raymar Oswald efficiency estimation for swept wings
def induced_drag(AR, LE_sweep, C_L, flap_deflection, h_winglet, b, alt): # b is the wing span of the jet, alt is the height of the aircraft above the ground, h_winglet is the height of the winglet
    delta_AR = 1.9 * (h_winglet / b) * AR
    AR = AR + delta_AR
    e = 4.6 * (1 - 0.045 * AR ** 0.68) * (math.cos(LE_sweep)) ** 0.15 - 3.1 
    delta_e = 0.0026 * flap_deflection
    e = e * delta_e
    K = 1 / (pi * e * AR)
    if alt < b/2:
        K = 33 * (alt/b) ** 1.5 / (1 + 33 * (alt/b) ** 1.5)
    c_induced = K * C_L ** 2
    return c_induced, e, AR

## wave drag
## min cp0 at 1.5 alpha -1.02
##M_crit = 0.648
def wave_C_D(Mach, M_DD):
    C_D_wave = 0.002*(1+2.5*(M_DD-Mach)/0.05)**(-1)
    return C_D_wave



def final_C_D_increase(C_D):
    C_D_final = C_D *1.03
    return C_D_final



