import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from XFLR import y_span, chord0, Ai, Cl0, ICd, Cm0
from XFLR import y_span10, chord10, Ai10, Cl10, ICd10, Cm10
from Integration import aoa_deg_abs
from scipy.optimize import curve_fit


from scipy import integrate, interpolate
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from shear_centre_location import shear_center_non_dim
import conditions as c
from Integration import x_grid


# --- Variables ---
#Wing
b_half = 9.7925 #half of the wing span in m
b = 19.585 #Span of the wing in m
C_r = 2.874 # Root cord in m
C_t = 1.04326 # Tip cord in m
#aoa_deg_abs = 0.0   # Angle of attack in degrees

#conditions
g = 9.81 # Gravitational constant m/s^2
V_inf = c.velocity  # Freestream velocity in m/s
rho   = c.density # Air density in kg/m^3
#Mass
M_wing = 932.9 # mass of the wing in kg
M_fuel_T1 = 533.6655 # mass of fuel in fuel tank 1 (close to the fuselage) in kg
M_fuel_T2 = 881.2825 # mass of fuel in fuel tank 2 (after landing gear) in kg
W_main_gear = 1245.87 #weight of landing gear (already accounted for there being two weight split half per wing (already halved)) in N



# --- Cord lengths ---------
def cordlength(tip_cord, root_cord, fraction_half_span):
    
    cord = fraction_half_span * (tip_cord - root_cord) + root_cord

    return cord

C_19 = cordlength(C_t, C_r, 0.19)
C_24 = cordlength(C_t, C_r, 0.24)
C_90 = cordlength(C_t, C_r, 0.9)



# --- determining loading functions ----------------------------------------------------------------




n = len(Cl0)
y_tip = y_span[n//2:]        # last half of spanwise locations
Cl0_tip = Cl0[n//2:] 
Cl10_tip = Cl10[n//2:] 
ICd0_tip = ICd[n//2:]  
ICd10_tip = ICd10[n//2:] 
Cm0_tip = Cm0[n//2:]  
Cm10_tip = Cm10[n//2:]          # corresponding Cl
chord_tip = chord0[n//2:]     # corresponding chord lengths

"""To produce NVM for AOA 10 degrees, uncomment THIS!
# --- FOR AOA 10!!!!---- 
L_prime = compute_lift_line_load(chord10, Cl10, V_inf, rho)
D_prime = compute_drag_line_load(chord10, ICd10, V_inf, rho)
N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg_abs=10.0)
M_prime = compute_section_moment_density(chord10, Cm10, V_inf, rho)
y, q_func, d_func, t_func = build_q_d_t_functions(y_span10, chord10, N_prime, M_prime, 10.43, 0.3, 0.7, 0.25)
#ALSO for AOA 10!!! 
n = len(Cl10)
y_tip = y_span10[n//2:]        # last half of spanwise locations
Cl_tip = Cl10[n//2:]           # corresponding Cl
chord_tip = chord10[n//2:]     # corresponding chord lengths
#"""

# --- Interpolating to 500 points ---
n_points = 500
y_interp = np.linspace(y_tip.min(), y_tip.max(), n_points)
Cl0_interp = np.interp(y_interp, y_tip, Cl0_tip)
Cl10_interp = np.interp(y_interp, y_tip, Cl10_tip)
ICd0_interp = np.interp(y_interp, y_tip, ICd0_tip)
ICd10_interp = np.interp(y_interp, y_tip, ICd10_tip)
Cm0_interp = np.interp(y_interp, y_tip, Cm0_tip)
Cm10_interp = np.interp(y_interp, y_tip, Cm10_tip)
chord_interp = np.interp(y_interp, y_tip, chord_tip)

# --- Lift line load function I just copied Caia's---
def compute_lift_line_load(chord: np.ndarray,
                           Cl0: np.ndarray,
                           Cl10: np.ndarray,
                           aoa_deg_abs: float, 
                           V_inf: float,
                           rho: float = 1.225) -> np.ndarray:
    """
    Computes spanwise lift per unit length:
        L'(y) = 0.5 * rho * V^2 * Cl(y) * c(y)
    """
    Cl = (Cl10 - Cl0) / 10 * aoa_deg_abs + Cl0

    q_inf = 0.5 * rho * V_inf**2
    return q_inf * chord * Cl

def compute_drag_line_load(chord: np.ndarray, ICd0: np.ndarray, ICd10: np.ndarray, aoa_deg_abs: float, V_inf: float, rho: float = 1.225) -> np.ndarray:
    """
    D'(y) = 0.5 * rho * V^2 * Cd(y) * c(y)
    """
    ICd = (ICd10 - ICd0) / 10 * aoa_deg_abs + ICd0

    q_inf = 0.5 * rho * V_inf**2
    D_prime = q_inf * chord * ICd
    return D_prime

L_prime = compute_lift_line_load(chord_interp, Cl0_interp, Cl10_interp, aoa_deg_abs, V_inf)
D_prime = compute_drag_line_load(chord_interp, ICd0_interp, ICd10_interp, aoa_deg_abs, V_inf)


def compute_normal_force_distribution(L_prime: np.ndarray,
                                      D_prime: np.ndarray,
                                      aoa_deg_abs: float) -> np.ndarray:
    """
    N'(y) = L'(y) * cos(aoa) + D'(y) * sin(aoa)
    """
    aoa_rad = np.radians(aoa_deg_abs)
    N_prime = L_prime * np.cos(aoa_rad) + D_prime * np.sin(aoa_rad)
    return N_prime

def compute_section_moment_density(chord: np.ndarray,
                                   Cm0: np.ndarray,
                                   Cm10: np.ndarray,
                                   aoa_deg_abs: float,
                                   V_inf: float,
                                   rho: float = 1.225) -> np.ndarray:
    """
    M'(y) = Cm(y) * q_inf * c(y)^2 
    """

    Cm = (Cm10 - Cm0) / 10 * aoa_deg_abs + Cm0

    q_inf = 0.5 * rho * V_inf**2
    M_prime = Cm * q_inf * chord**2
    return M_prime




# --- Compute L'(y) for the tip half ---
  # freestream velocity [m/s]
N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg_abs)
M_prime = compute_section_moment_density(chord_interp, Cm0_interp, Cm10_interp, aoa_deg_abs, V_inf, rho)



"""# --- Plotting the lift distribution to check if sensical---
plt.plot(y_interp, L_prime)
plt.xlabel("Spanwise position y [m]")
plt.ylabel("Lift per unit span L' [N/m]")
plt.title("Spanwise Lift Distribution (Tip Half)")
plt.grid(True)
plt.show()
"""

# --- Create grid ---
y_vals = x_grid
N = len(y_vals) # we choose this to be 500

# --- Index boundaries ---
i_19 = int(0.19 * N)   # Tank 1:   0% → 19%, what this does it returns an index position for what 19% would be.
i_24 = int(0.24 * N)   # Tank 2:  24% → 90%
i_90 = int(0.90 * N) 


# --- Wing structure load ---
def wing_weight_distribution(mass_wing, grav_const, wing_span, tip_cord, root_cord, y_values, aoa_deg_abs):
    aoa_rad = np.radians(aoa_deg_abs)
    return np.cos(aoa_rad)*(2 * mass_wing * grav_const) / (wing_span * (tip_cord + root_cord)) * \
           (root_cord + (2 * (tip_cord - root_cord) / wing_span) * y_values)

# --- Fuel tanks ---
def Fuel_distribution_tank_1(mass_fuel, grav_const, wing_span, root_cord, cord_19, y_values, aoa_deg_abs):
    aoa_rad = np.radians(aoa_deg_abs)
    return np.cos(aoa_rad)*(4 * mass_fuel * grav_const) / (0.19 * wing_span * (root_cord + cord_19)) * \
           (root_cord + ((cord_19 - root_cord) / (0.095 * wing_span)) * y_values)

def Fuel_distribution_tank_2(mass_fuel, grav_const, wing_span, cord_24, cord_90, y_values, aoa_deg_abs):
    aoa_rad = np.radians(aoa_deg_abs)
    return np.cos(aoa_rad)*(4 * mass_fuel * grav_const) / (0.66 * wing_span * (cord_24 + cord_90)) * \
           (cord_24 + ((cord_90 - cord_24) / (0.33 * wing_span)) * y_values)



# --- Combined load array ---
combined_loads = np.zeros_like(y_vals)

# --- Structural load over full span ---

wing_weight_only = wing_weight_distribution(M_wing, g, b, C_t, C_r, y_vals, aoa_deg_abs)


# --- Tank 1 load (0% → 19%) ---
y_t1 = y_vals[:i_19]                      # local y inside tank 1 region
W_t1 = Fuel_distribution_tank_1(M_fuel_T1, g, b, C_r, C_19, y_t1, aoa_deg_abs)

# --- Tank 2 load (24% → 90%) ---
# Shift y so Tank 2 formula sees y=0 at 24%
y_t2_local = y_vals[i_24:i_90] - y_vals[i_24]
W_t2 = Fuel_distribution_tank_2(M_fuel_T2, g, b, C_24, C_90, y_t2_local, aoa_deg_abs)



# Spanwise segment for the gear
y_gear = y_vals[107]                  # local y inside gear region
gear_load_per_point = W_main_gear / 0.097925 # distribut   ed the landing gear weight over 10 cm



""" --- Adding the distributions --------------------------------------------------------------------------
#----------------If AoA is 0 deg ------------------
combined_loads[:] -= wing_weight_only                                        # struc Aoa=0
combined_loads[:i_19] -= W_t1                                                # Tank 1 AoA=0
combined_loads[i_24:i_90] -= W_t2                                            # Tank 2 AoA=0
combined_loads[105:109] -= gear_load_per_point                             #Landing gear AoA=0
combined_loads[:] += N_prime                                                 #Lift AoA=0
#"""


#"""----------------Combined loads ------------------           # to use this add a # in front of this line
combined_loads[:] -= wing_weight_only * np.cos(np.deg2rad(aoa_deg_abs))               # struc
combined_loads[:i_19] -= W_t1 * np.cos(np.deg2rad(aoa_deg_abs))                       # Tank 1
combined_loads[i_24:i_90] -= W_t2 * np.cos(np.deg2rad(aoa_deg_abs))                   # Tank 2
combined_loads[i_19:i_24] -= gear_load_per_point * np.cos(np.deg2rad(aoa_deg_abs))    #Landing gear
combined_loads[:] += N_prime                                                      #Lift & Drag
#combined_loads *= abs(c.load_factor)
#"""

# --- SHEAR FORCE S(y) -----------------------------------------------------------------------------------
# Integrate q(y) from tip -> root
S_vals_tip_to_root = cumulative_trapezoid(combined_loads[::-1], y_vals[::-1], initial=0)

# Flip so that y increases from root → tip
S_vals = S_vals_tip_to_root[::-1]

# --- BENDING MOMENT M(y) ----------------------------------------------------------------------------------
# Integrate S(y) from tip -> root
M_vals_tip_to_root = cumulative_trapezoid(S_vals[::-1], y_vals[::-1], initial=0)

# Flip back
M_vals = M_vals_tip_to_root[::-1]





# ---------- PLOTS ---------------------------------------------------------------------------------------------

plt.figure(figsize=(10,10))

# q(y)
plt.subplot(3,1,1)
plt.plot(y_vals, combined_loads)
plt.xlabel("Spanwise Position y [m]")
plt.ylabel("q(y) [N/m]")
plt.title("Distributed Load q(y)")
plt.grid(True)

# S(y)+
plt.subplot(3,1,2)
plt.plot(y_vals, S_vals)
plt.xlabel("Spanwise Position y [m]")
plt.ylabel("Shear Force S(y) [N]")
plt.title("Shear Force Diagram")
plt.grid(True)

# M(y)
plt.subplot(3,1,3)
plt.plot(y_vals, M_vals)
plt.xlabel("Spanwise Position y [m]")
plt.ylabel("Bending Moment M(y) [Nm]")
plt.title("Bending Moment Diagram")
plt.grid(True)

plt.tight_layout()
plt.show()

