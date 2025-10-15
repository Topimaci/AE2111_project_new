

#Values we choose from assumptions:
#import * from Planform_DESIGN1 

#-------------------Import Values from Iterative procedure

S_w = 49.63
b_w = 22.278
c_w = 2.39
c_r = 3.27
MTOM = 12520
fuel_mass_fraction = 0.333
#c_wn is chord mac of main wing, c_r is root chord

#--------------------Assumptions and Constants needed for the Horizontal and Vertical tail surface area calculation-----------------------------------#

# Aspect Ratio horizontal tail, ranges from 3 to 4 recommended*

AR_h = 3.5 

# Taper Ratio horizontal tail, ranges from 0.6 to 1 for T-tail. 

taper_h = 0.7

#10.4 to 50, where 10.4 originates from final planform WP2, being the lower limit (adsee book)

sweep_leading_edge_v = 20

# Coefficient of volume, also chosen based on aircraft type: for business jet, 0.059 to 0.093

V_v = 0.06


#nacelle length excludes the cone
fus_length = 15.2
nacelle_length = 0.48


#-------------------------CG_aft position calculation in order to determine arm for Vertical and Horizontal tail area-----------------------#

#  Mass Fraction Estimations 

m_OE = 0.607
m_wing = 0.135
m_fus = 0.105
m_t = 0.043
m_eng = 0.1
m_nac = 0.018
m_lg = 0.036
m_fe = 0.17
m_unacc = 0


#sanity check available to check the mass fractions equate to 0 total. 

m_total = m_OE - (m_wing + m_fus + m_t + m_eng + m_nac + m_lg + m_fe + m_unacc)

print("Mass fractions total (should be equal to 0 if all is correct:" , m_total)

#-------------------------Fuselage and Wing Assembly moment arm calculations------------------------------------#


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


#Individual component moments

mom_fus = arm_fus * m_fus
mom_t = arm_t * m_t
mom_eng = arm_eng * m_eng
mom_nacelle = arm_nac * m_nac
mom_fe = arm_fe * m_fe

mom_wing = arm_wing * m_wing


#Total moment fuselage group

mom_fuselage_group = mom_fus + mom_t + mom_eng + mom_nacelle + mom_fe


#Total moment wing group

mom_wing_group = mom_wing




#constants

xc_OEWCG = 0.25 #assumption
cg_payload = 8 #assumed that it is located in the middle of the cabin (l_cab = 10m) + length_nose = 3m
m_payload = 1010 #kg

#functions

def Xlemac(m_fuselage_group, mom_fuselage_group, mom_wing_group, xc_OEWCG, c_w):
    Xlemac = mom_fuselage_group / m_fuselage_group + c_w * ((mom_wing_group/m_fuselage_group) * 0.4 - xc_OEWCG * (1 + mom_wing_group/m_fuselage_group)) #0.4 is an assumption and corresponds to (x/c_mac)_WACG, essentially location of the wing assembly center of gravity with respect to the mean aerodynamic chord
    return Xlemac

def X_OEM(Xlemac, c_w, xc_OEWCG):
    X_OEM = Xlemac + c_w * xc_OEWCG
    return X_OEM

def cgaft(MTOM, X_OEM, fuel_mass_fraction, Xlemac):

    CGofOEMandMAXPAYLOAD = ((MTOM * m_OE) * X_OEM + cg_payload * m_payload) / ((MTOM * m_OE) + m_payload)
    
    CGofOEMandMAXPAYLOADandFUEL = ((MTOM * m_OE) * X_OEM + (cg_payload * m_payload) + (fuel_mass_fraction * MTOM) * (Xlemac + c_w * 0.15)) / ((MTOM * m_OE) + m_payload + fuel_mass_fraction * MTOM) #0.15 is an assumption of center of gravity of the fuel with respect to the mac
    
    CGofOEMandFUEL = ((MTOM * m_OE) * X_OEM + (Xlemac + c_w * 0.15) * (fuel_mass_fraction * MTOM)) / ((MTOM * m_OE) + (fuel_mass_fraction * MTOM))
    
    cgaft = max(CGofOEMandMAXPAYLOAD, CGofOEMandMAXPAYLOADandFUEL, CGofOEMandFUEL)
    

    print("1", CGofOEMandMAXPAYLOAD)
    print("2", CGofOEMandMAXPAYLOADandFUEL)
    print("3", CGofOEMandFUEL)


    return cgaft




def calculate_surface_area_vertical_tail(V_v, S_w, b_w, l_v):

    S_v = ( V_v * S_w * b_w ) / l_v

    return S_v



V_h = 0.61
#coefficient of volume, also chosen based on aircraft type: for business jet, 0.51 to 0.99


def calculate_surface_area_horizontal_tail(V_h, S_w, c_w, l_h):

    S_h = ( V_h * S_w * c_w ) / l_h

    return S_h


Xlemac = Xlemac(m_fuselage_group, mom_fuselage_group, mom_wing_group, xc_OEWCG, c_w)
X_OEM = X_OEM(Xlemac, c_w, xc_OEWCG)
cgaft = cgaft(MTOM, X_OEM, fuel_mass_fraction, Xlemac)
print("cg", cgaft)
l_v = 0.9*fus_length - cgaft
l_h = l_v 

print(l_v)

S_v = calculate_surface_area_vertical_tail(V_v, S_w, b_w, l_v)
S_h = calculate_surface_area_horizontal_tail(V_h, S_w, c_w, l_h)



print("Vertical wing area: ", S_v)
print("Horizontal wing area: ", S_h)


