import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from XFLR import y_span, chord, Ai, Cl, ICd, Cm


from scipy import integrate, interpolate
from scipy.integrate import cumulative_trapezoid



from TorqueDist import compute_lift_line_load
from TorqueDist import build_q_d_t_functions
from TorqueDist import compute_normal_force_distribution
from TorqueDist import compute_drag_line_load
from TorqueDist import compute_section_moment_density

# --- Variables ---
#Wing
b_half = 9.7925 #half of the wing span in m
b = 19.585 #Span of the wing in m
C_r = 2.874 # Root cord in m
C_t = 1.04326 # Tip cord in m
aoa_deg = 0.0   # Angle of attack in degrees

#conditions
g = 9.81 # Gravitational constant m/s^2
V_inf = 10.0  # Freestream velocity in m/s
rho   = 1.225 # Air density in kg/m^3

#Mass
M_wing = 932.9 # mass of the wing in kg
M_fuel_T1 = 533.6655 # mass of fuel in fuel tank 1 (close to the fuselage) in kg
M_fuel_T2 = 881.2825 # mass of fuel in fuel tank 2 (after landing gear) in kg


# --- Cord lengths ---------
def cordlength(tip_cord, root_cord, fraction_half_span)
    
    cord = fraction_half_span * (tip_cord - root_cord) + root_cord

    return cord

C_19 = cordlength(C_t, C_r, 0.19)
C_24 = cordlength(C_t, C_r, 0.24)
C_90 = cordlength(C_t, C_r, 0.9)


# --- determining loading functions ----------------------------------------------------------------
# --- Lift ---
L_prime = compute_lift_line_load(chord, Cl, V_inf, rho)
D_prime = compute_drag_line_load(chord, ICd, V_inf, rho)
N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg)
M_prime = compute_section_moment_density(chord, Cm, V_inf, rho)
y, q_func, d_func, t_func = build_q_d_t_functions(y_span, N_prime, M_prime)

# Create y-grid only over the half-span
y_vals = np.linspace(0, b_half, 500)

# Evaluate q(y)
q_vals = q_func(y_vals)

# --- Wing weight --- 

#Wing weight distribution function
def wing_weight_distribution(mass_wing, grav_const, wing_span, tip_cord, root_cord, y_values):
    
    W_struc_func = (2 * mass_wing * grav_const)/(wing_span(tip_cord + root_cord)) * (root_cord + ((2*(tip_cord - root_cord))/wing_span) * y_values)


    return W_struc_func

# --- Fuel weight: Tank 1 and Tank 2 ---------------------

def Fuel_distribution_tank_1(mass_fuel, grav_const, wing_span, root_cord, cord_19, y_values):
    
    W_fuel_tank_1_func = ((4*mass_fuel*grav_const)/(0.19 * wing_span * (root_cord + cord_19))) * (root_cord + (((cord_19 - root_cord)/0.19*wing_span) * y_values)) 
    # Function from 0% of the half wingspan to 19% of the half wingspan, not to be used outside of this range
    return W_fuel_tank_1_func

def Fuel_distribution_tank_2(mass_fuel, grav_const, wing_span, cord_24, cord_90, y_values):
    
    W_fuel_tank_2_func = ((4*mass_fuel*grav_const)/(0.66 * wing_span * (cord_24 + cord_90))) * (cord_24 + (((cord_90 - cord_24)/0.66*wing_span) * y_values)) 
    # Function from 24% of the half wingspan to 90% of the half wingspan, not to be used outside of this range
    return W_fuel_tank_2_func


import numpy as np

def wing_weight_distribution(...):
    ...
    return W_struc_func

def Fuel_distribution_tank_1(...):
    ...
    return W_fuel_tank_1_func

def Fuel_distribution_tank_2(...):
    ...
    return W_fuel_tank_2_func


# y_values: array from 0 → half-span
y = y_values

# Empty array for combined loads
combined_loads = np.zeros_like(y)

# Structural load (full span)
combined_loads += wing_weight_distribution(...)

# Tank 1: only 0%–19% of half-span
mask_t1 = (y >= 0) & (y <= 0.19 * wing_span/2)
combined_loads[mask_t1] += Fuel_distribution_tank_1(...)[mask_t1]

# Tank 2: only 24%–90% of half-span
mask_t2 = (y >= 0.24 * wing_span/2) & (y <= 0.90 * wing_span/2)
combined_loads[mask_t2] += Fuel_distribution_tank_2(...)[mask_t2]


# --- SHEAR FORCE S(y) -----------------------------------------------------------------------------------
# Integrate q(y) from tip -> root
S_vals_tip_to_root = cumulative_trapezoid(q_vals[::-1], y_vals[::-1], initial=0)

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
plt.plot(y_vals, q_vals)
plt.xlabel("Spanwise position y [m]")
plt.ylabel("q(y) [N/m]")
plt.title("Distributed Load q(y)")
plt.grid(True)

# S(y)
plt.subplot(3,1,2)
plt.plot(y_vals, S_vals)
plt.xlabel("Spanwise position y [m]")
plt.ylabel("Shear Force S(y) [N]")
plt.title("Shear Force Diagram")
plt.grid(True)

# M(y)
plt.subplot(3,1,3)
plt.plot(y_vals, M_vals)
plt.xlabel("Spanwise position y [m]")
plt.ylabel("Bending Moment M(y) [Nm]")
plt.title("Bending Moment Diagram")
plt.grid(True)

plt.tight_layout()
plt.show()

