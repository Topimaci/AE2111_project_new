import numpy as np

import matplotlib.pyplot as plt


#---Constants---------------------------------------------------------------

# --- Design 1 ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_S12_D1 = 12.0
Stringers_S34_D1 = 8.0

A_D1 = 0.0002 # Area stringer (crossection) in m^2

t_skin_D1 = 0.005 # Thickness of the skin in m

# --- Design 2 ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_S12_D2 = 13.0
Stringers_S34_D2 = 8.0

A_D2 = 0.00028 # Area stringer (crossection) in m^2

t_skin_D2 = 0.003 # Thickness of the skin in m

# --- Design 3 ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_S12_D3 = 8.0
Stringers_S34_D3 = 4.0

A_D3 = 0.00022 # Area stringer (crossection) in m^2

t_skin_D3 = 0.005 # Thickness of the skin in m

# --- Other constants ---
E = 71 * (10**9) # Young's modulus in Pa
k_S123 = 4 # Boundary condition stringer section 12 constant (=4 due to the assumption that both ends are clamped)
k_S4 = 0.25 # boundary condition stringer section 34(fixed-free)
wing_half = 9.79 # half of the wingspan in m
L = 0.25 * wing_half # Length of a  uninterupted stringer = half of half the wingspan since # of stringers change at half the wingspan 


# --- Geometry of stringer (L-shaped) -------------------------------------------
# Assumptions: t = t_skin, h = 2w, double thickness contribution in intersection is neglegted(due to thin-walled)
def Stringer_Geometry(Area_stringer, t_skin): # area in m^2, t in m
    
    width = (Area_stringer + t_skin**2) / ( 2 * t_skin)
    height = width

    return width, height # in m

w_stinger_D1, h_stringer_D1 = Stringer_Geometry(A_D1, t_skin_D1)
w_stinger_D2, h_stringer_D2 = Stringer_Geometry(A_D2, t_skin_D2)
w_stinger_D3, h_stringer_D3 = Stringer_Geometry(A_D3, t_skin_D3)
print(w_stinger_D1, h_stringer_D1)

# --- Moment of Inertia of stringer ----------------------------
#Assumptions: thinwalled approx 
def MoI_stringer(b, h, t):
    y_bar = (h**2 + b*t - t**2) / (2*(b + h - t))

    I_stringer = (t*h**3)/12 + t*h*(y_bar - h/2)**2 + (b*t**3)/12 + b*t*(y_bar - t/2)**2 - t**4/12 - t**2*(y_bar - t/2)**2

    return I_stringer # in m^4
    
I_stringer_D1 = MoI_stringer(w_stinger_D1, h_stringer_D1, t_skin_D1)
I_stringer_D2 = MoI_stringer(w_stinger_D2, h_stringer_D2, t_skin_D2)
I_stringer_D3 = MoI_stringer(w_stinger_D3, h_stringer_D3, t_skin_D3)
#print(I_stringer_D1)

# ---Critical stress--------------------------------------------------------
# The critical stress for a certain stringer

def column_critical_stress(k, E, I, length, Area):

    sigma_crit = (k * (np.pi)**2 * E * I )/ ((length**2) * Area )

    return sigma_crit

sigma_crit_per_stringer_S1_D1 = column_critical_stress(k_S123, E, I_stringer_D1, L, A_D1)
sigma_crit_per_stringer_S2_D1 = column_critical_stress(k_S123, E, I_stringer_D1, L, A_D1)
sigma_crit_per_stringer_S3_D1 = column_critical_stress(k_S123, E, I_stringer_D1, L, A_D1)
sigma_crit_per_stringer_S4_D1 = column_critical_stress(k_S4, E, I_stringer_D1, L, A_D1)

sigma_crit_per_stringer_S1_D2 = column_critical_stress(k_S123, E, I_stringer_D2, L, A_D2)
sigma_crit_per_stringer_S2_D2 = column_critical_stress(k_S123, E, I_stringer_D2, L, A_D2)
sigma_crit_per_stringer_S3_D2 = column_critical_stress(k_S123, E, I_stringer_D2, L, A_D2)
sigma_crit_per_stringer_S4_D2 = column_critical_stress(k_S4, E, I_stringer_D2, L, A_D2)

sigma_crit_per_stringer_S1_D3 = column_critical_stress(k_S123, E, I_stringer_D3, L, A_D3)
sigma_crit_per_stringer_S2_D3 = column_critical_stress(k_S123, E, I_stringer_D3, L, A_D3)
sigma_crit_per_stringer_S3_D3 = column_critical_stress(k_S123, E, I_stringer_D3, L, A_D3)
sigma_crit_per_stringer_S4_D3 = column_critical_stress(k_S4, E, I_stringer_D3, L, A_D3)
#print(sigma_crit_per_stringer_S12)
#print(sigma_crit_per_stringer_S34)


# in array form * the # of stringers
sigma_crit_column_S1_D1 = np.full(125, sigma_crit_per_stringer_S1_D1) * Stringers_S12_D1
sigma_crit_column_S2_D1 = np.full(125, sigma_crit_per_stringer_S2_D1) * Stringers_S12_D1
sigma_crit_column_S3_D1 = np.full(125, sigma_crit_per_stringer_S3_D1) * Stringers_S34_D1
sigma_crit_column_S4_D1 = np.full(125, sigma_crit_per_stringer_S4_D1) * Stringers_S34_D1
sigma_crit_column_S12_D1 = np.append(sigma_crit_column_S1_D1,sigma_crit_column_S2_D1) 
sigma_crit_column_S34_D1 = np.append(sigma_crit_column_S3_D1, sigma_crit_column_S4_D1)
sigma_crit_column_full_D1 = np.append(sigma_crit_column_S12_D1,sigma_crit_column_S34_D1)

sigma_crit_column_S1_D2 = np.full(125, sigma_crit_per_stringer_S1_D2) * Stringers_S12_D2
sigma_crit_column_S2_D2 = np.full(125, sigma_crit_per_stringer_S2_D2) * Stringers_S12_D2
sigma_crit_column_S3_D2 = np.full(125, sigma_crit_per_stringer_S3_D2) * Stringers_S34_D2
sigma_crit_column_S4_D2 = np.full(125, sigma_crit_per_stringer_S4_D2) * Stringers_S34_D2
sigma_crit_column_S12_D2 = np.append(sigma_crit_column_S1_D2,sigma_crit_column_S2_D2)
sigma_crit_column_S34_D2 = np.append(sigma_crit_column_S3_D2,sigma_crit_column_S4_D2)
sigma_crit_column_full_D2 = np.append(sigma_crit_column_S12_D2,sigma_crit_column_S34_D2)


sigma_crit_column_S1_D3 = np.full(125, sigma_crit_per_stringer_S1_D3) * Stringers_S12_D3
sigma_crit_column_S2_D3 = np.full(125, sigma_crit_per_stringer_S2_D3) * Stringers_S12_D3
sigma_crit_column_S3_D3 = np.full(125, sigma_crit_per_stringer_S3_D3) * Stringers_S34_D3
sigma_crit_column_S4_D3 = np.full(125, sigma_crit_per_stringer_S4_D3) * Stringers_S34_D3
sigma_crit_column_S12_D3 = np.append(sigma_crit_column_S1_D3,sigma_crit_column_S2_D3)
sigma_crit_column_S34_D3 = np.append(sigma_crit_column_S3_D3,sigma_crit_column_S4_D3)
sigma_crit_column_full_D3 = np.append(sigma_crit_column_S12_D3,sigma_crit_column_S34_D3)


#'''
print(sigma_crit_column_full_D1)
print(sigma_crit_column_full_D2)
print(sigma_crit_column_full_D3)
#'''

#'''
x = np.linspace(0, 10, 500)
plt.plot(x, sigma_crit_column_full_D1)
plt.title("D1")
plt.xlabel("x axis")
plt.ylabel("crit sigma")
plt.show()

plt.plot(x, sigma_crit_column_full_D2)
plt.title("D2")
plt.xlabel("x axis")
plt.ylabel("crit sigma")
plt.show()

plt.plot(x, sigma_crit_column_full_D3)
plt.title("D3")
plt.xlabel("x axis")
plt.ylabel("crit sigma")
plt.show()
#'''