import numpy as np


#---Constants---------------------------------------------------------------

# --- Design ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_S1 = 12.0
Stringers_S2 = 12.0
Stringers_S3 = 8.0
Stringers_S4 = 8.0

A = 0.0002 # Area stringer (crossection) in m^2

t_skin = 0.005 # Thickness of the skin in m

# --- Other constants ---
E = 71 * (10**9) # Young's modulus in Pa
k = 4 # Boundary condition constant (=4 due to the assumption that both ends are clamped)
wing_half = 9.79 # half of the wingspan in m
L = 0.5 * wing_half # Length of a  uninterupted stringer = half of half the wingspan since # of stringers change at half the wingspan 


# --- Geometry of stringer (L-shaped) -------------------------------------------
# Assumptions: t = t_skin, h = 2w, double thickness contribution in intersection is neglegted(due to thin-walled)
def Stringer_Geometry(Area_stringer, t_skin): # area in m^2, t in m
    
    width = (Area_stringer) / ( 3 * t_skin)
    height = width / 2

    return width, height # in m

w_stinger, h_stringer = Stringer_Geometry(A, t_skin)

# --- Moment of Inertia of stringer ----------------------------
#Assumptions: thinwalled approx 
def MoI_stringer(width, height, thickness):
    y_bar = (width + 0.5*(height**2))/ (width * height * thickness)

    I_stringer = (width * thickness * y_bar) + (thickness * (height**3))/12 + (height * thickness) * ((0.5*height - y_bar)**2)

    return I_stringer # in m^4
    
I_stringer = MoI_stringer(w_stinger, h_stringer, t_skin)


# ---Critical stress--------------------------------------------------------
# The critical stress for a certain stringer

def column_critical_stress(k, E, I, length, Area):

    sigma_crit = (k * (np.pi)**2 * E * I )/ ((length**2) * Area )

    return sigma_crit

sigma_crit_per_stringer = column_critical_stress(k, E, I_stringer, L, A)

# in array form
sigma_crit_column = np.full(500, sigma_crit_per_stringer)
print(sigma_crit_per_stringer)