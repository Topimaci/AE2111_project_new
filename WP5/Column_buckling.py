import numpy as np


#---Constants---------------------------------------------------------------

# --- Design ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_S1 = 6.0
Stringers_S2 = 6.0
Stringers_S3 = 3.0
Stringers_S4 = 3.0

A = 45.0 # Area stringer (crossection) in cm^2

t_skin = 5.0 # Thickness of the skin in mm

# --- Other constants ---
E = 71 * 10**9 # Young's modulus in Pa
k = 5.57 # Boundary condition constant (=4 due to the assumption that both ends are clamped)
wing_half = 9.79 # half of the wingspan in m


# --- Geometry of stringer (L-shaped) -------------------------------------------
# Assumptions: t = t_skin, h = 2w, double thickness contribution in intersection is neglegted(due to thin-walled)
def Stringer_Geometry(Area_stringer, t_skin): # area in cm^2, t in mm
    
    width = (Area_stringer * 10**2) / ( 3 * t_skin)
    height = width / 2

    return width, height # in mm

w_stinger, h_stringer = Stringer_Geometry(A, t_skin)

# --- Moment of Inertia of stringer ----------------------------
#Assumptions: thinwalled approx 
def MoI_stringer(width, height, thickness):
    y_bar = (width + 0.5*(height**2))/ (width * height * thickness)

    I_stringer = (width * thickness * y_bar) + (thickness * (height**3))/12 + (height * thickness) * ((0.5*height - y_bar)**2)

    return I_stringer # in mm^4
    
I_stringer = MoI_stringer(w_stinger, h_stringer, t_skin)


# --- Length of stringer ---------------------









# ---Critical stress---
# The critical stress for a certain stringer

def column_critical_stress(k, E, I, length, Area):

    sigma_crit = (k * (np.pi)**2 * E * (I * 10**-12 ) )/ ((length**2) * (Area * 10**-4))

    return sigma_crit
