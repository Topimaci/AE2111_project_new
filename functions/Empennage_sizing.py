def calculate_tail_surface_areas(
    S_w,
    b_w,
    c_w,
    c_r,
    MTOM,
    fuel_mass_fraction,
    #-------------------------Mass Fraction Estimations -------------------------#
    m_OE,
    m_wing,
    m_fus,
    m_t,
    m_eng,
    m_nac,
    m_fe,
):

     #-----------CONSTANTS--------#
    #Following values can be changed according to needs, but do not follow iterative procedure 
    fus_length = 15.2
    nacelle_length = 3
    # Coefficient of volume, also chosen based on aircraft type: for business jet
    V_v = 0.06
    V_h = 0.61
    xc_OEWCG = 0.25 # assumption
    cg_payload = 8  # assumed that it is located in the middle of the cabin
    m_payload = 1010 # kg

    # Total Mass fraction fuselage group
    m_fuselage_group = m_fus + m_t + m_eng + m_nac + m_fe 

    # Total Mass fraction wing group 
    m_wing_group = m_wing

    # Moment arm calculations for each subcomponent of Fuselage and Wing subgroups
    arm_fus = 0.4 * fus_length
    arm_t = 0.9 * fus_length
    arm_eng = 0.4 * nacelle_length  + 0.75 * fus_length
    arm_nac = 0.4 * nacelle_length  + 0.75 * fus_length
    arm_fe = 0.4 * fus_length
    arm_wing = 0.4 * c_r

    # Moment calculations for each subcomponent of Fuselage and Wing subgroups
    mom_fuselage_group = arm_fus * m_fus + arm_t * m_t + arm_eng * m_eng + arm_nac * m_nac + arm_fe * m_fe
    mom_wing_group = arm_wing * m_wing_group

    # functions
    def Xlemac_func(m_fuselage_group, mom_fuselage_group, mom_wing_group, xc_OEWCG, c_w):
        Xlemac = mom_fuselage_group / m_fuselage_group + c_w * (
            (mom_wing_group/m_fuselage_group) * 0.4 - 
            xc_OEWCG * (1 + mom_wing_group/m_fuselage_group)
        ) #0.4 is an assumption and corresponds to (x/c_mac)_WACG
        return Xlemac

    def X_OEM_func(Xlemac, c_w, xc_OEWCG):
        X_OEM = Xlemac + c_w * xc_OEWCG
        return X_OEM

    def cgaft_func(MTOM, X_OEM, fuel_mass_fraction, Xlemac):
        CGofOEMandMAXPAYLOAD = ((MTOM * m_OE) * X_OEM + cg_payload * m_payload) / ((MTOM * m_OE) + m_payload)
        CGofOEMandMAXPAYLOADandFUEL = ((MTOM * m_OE) * X_OEM + (cg_payload * m_payload) + 
                                       (fuel_mass_fraction * MTOM) * (Xlemac + c_w * 0.15)) / ((MTOM * m_OE) + m_payload + fuel_mass_fraction * MTOM) #0.15 assumption
        CGofOEMandFUEL = ((MTOM * m_OE) * X_OEM + (Xlemac + c_w * 0.15) * (fuel_mass_fraction * MTOM)) / ((MTOM * m_OE) + (fuel_mass_fraction * MTOM))
        return max(CGofOEMandMAXPAYLOAD, CGofOEMandMAXPAYLOADandFUEL, CGofOEMandFUEL)

    def calculate_surface_area_vertical_tail(V_v, S_w, b_w, l_v):
        S_v = ( V_v * S_w * b_w ) / l_v
        return S_v

    def calculate_surface_area_horizontal_tail(V_h, S_w, c_w, l_h):
        S_h = ( V_h * S_w * c_w ) / l_h
        return S_h

    Xlemac = Xlemac_func(m_fuselage_group, mom_fuselage_group, mom_wing_group, xc_OEWCG, c_w)
    X_OEM = X_OEM_func(Xlemac, c_w, xc_OEWCG)
    cgaft_val = cgaft_func(MTOM, X_OEM, fuel_mass_fraction, Xlemac)

    # Tail arm distances
    l_v = 0.9*fus_length - cgaft_val
    l_h = l_v 

    # Tail surface areas
    S_v = calculate_surface_area_vertical_tail(V_v, S_w, b_w, l_v)
    S_h = calculate_surface_area_horizontal_tail(V_h, S_w, c_w, l_h)

    return S_v, S_h, l_v


# Example usage:
S_v, S_h = calculate_tail_surface_areas()
print("Vertical tail area:", S_v)
print("Horizontal tail area:", S_h)

S_v, S_h = calculate_tail_surface_areas(S_w, b_w, c_w, c_r, MTOM, fuel_mass_fraction, m_OE, m_wing, m_fus, m_t, m_eng, m_nac, m_lg, m_fe, m_unacc)


# Aspect Ratio horizontal tail, ranges from 3 to 4 recommended*
AR_h = 3.5
# Taper Ratio horizontal tail, ranges from 0.6 to 1 for T-tail.
taper_h = 0.7
#10.4 to 50, where 10.4 originates from final planform WP2
sweep_leading_edge_v = 20