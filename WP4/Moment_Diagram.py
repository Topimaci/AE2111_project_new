import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from XFLR import y_span, chord, Ai, Cl, ICd, Cm


from scipy import integrate, interpolate


from TorqueDist import compute_lift_line_load
from TorqueDist import build_q_d_t_functions
from TorqueDist import compute_normal_force_distribution
from TorqueDist import compute_drag_line_load
from TorqueDist import compute_section_moment_density

#Variables
b_half = 9.7925 #half of the wing span in m
V_inf = 10.0  # Freestream velocity in m/s
rho   = 1.225 # Air density in kg/m^3
aoa_deg = 0.0   # Angle of attack in degrees


#determining q_func
L_prime = compute_lift_line_load(chord, Cl, V_inf, rho)
D_prime = compute_drag_line_load(chord, ICd, V_inf, rho)
N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg)
M_prime = compute_section_moment_density(chord, Cm, V_inf, rho)
y, q_func, d_func, t_func = build_q_d_t_functions(y_span, N_prime, M_prime)


# Create y-grid only over the half-span
y_vals = np.linspace(0, b_half, 500)

# Evaluate q(y)
q_vals = q_func(y_vals)

# Plot
plt.figure(figsize=(8,4))
plt.plot(y_vals, q_vals)
plt.xlabel("Spanwise position y [m]")
plt.ylabel("q(y)")
plt.title("q(y) Over Half Span")
plt.grid(True)
plt.show()





# Integrating load line -> shear load line S(y)

#S, error_S = sp.integrate.quad(q_func, 0, b_half)

#print(S)

# Integrating S(y) -> M(y)      negative sign???
#M, error_M = sp.integrate.quad(S, 0, b_half)

#print(S)