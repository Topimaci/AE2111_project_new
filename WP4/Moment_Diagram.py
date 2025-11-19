import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from XFLR import y_span, chord, Ai, Cl, ICd, Cm


from scipy import integrate, interpolate


from TorqueDist import compute_lift_line_load
from TorqueDist import build_q_d_t_functions
from TorqueDist import N_prime
from TorqueDist import M_prime

#Variables
b_half = 9.7925 #half of the wing span in m
V_inf = 10.0  # Freestream velocity in m/s
rho   = 1.225 # Air density in kg/m^3

#determining q_func

y, q_func, d_func, t_func = build_q_d_t_functions(y_span, N_prime, M_prime)

# Integrating load line -> shear load line S(y)

#S, error_S = sp.integrate.quad(q_func, 0, b_half)

# Integrating S(y) -> M(y)      negative sign???
#M, error_M = sp.integrate.quad(S, 0, b_half)

#print(S)