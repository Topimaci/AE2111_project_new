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



# --- determining loading functions ---
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




# --- SHEAR FORCE S(y) --- 
# Integrate q(y) from tip -> root
S_vals_tip_to_root = cumulative_trapezoid(q_vals[::-1], y_vals[::-1], initial=0)

# Flip so that y increases from root → tip
S_vals = S_vals_tip_to_root[::-1]

# --- BENDING MOMENT M(y) ---
# Integrate S(y) from tip -> root
M_vals_tip_to_root = cumulative_trapezoid(S_vals[::-1], y_vals[::-1], initial=0)

# Flip back
M_vals = M_vals_tip_to_root[::-1]


# ---------- PLOTS ----------

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

