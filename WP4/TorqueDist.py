import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


# Variables

V_inf = 10.0  # Freestream velocity in m/s
rho   = 1.225 # Air density in kg/m^3
aoa = 0.0   # Angle of attack in degrees

# Importing data and computing aerodynamic line load

filename = 'angleofattack0.txt'

# 1. loading data

def load_xflr_data(filename: str):
    data = np.loadtxt(filename, skiprows=21)

    y_span = data[:, 0]
    chord  = data[:, 1]
    Ai     = data[:, 2]
    Cl     = data[:, 3]
    ICd    = data[:, 5]
    Cm     = data[:, 7]

    return y_span, chord, Ai, Cl, ICd, Cm


def compute_lift_line_load(chord: np.ndarray,
                           Cl: np.ndarray,
                           V_inf: float = 10.0,
                           rho: float = 1.225) -> np.ndarray:
    """
    L'(y) = 0.5 * rho * V^2 * Cl(y) * c(y)
    """
    q_inf = 0.5 * rho * V_inf**2
    L_prime = q_inf * chord * Cl
    return L_prime

def compute_drag_line_load(chord: np.ndarray,
                           ICd: np.ndarray,
                           V_inf: float = 10.0,
                           rho: float = 1.225) -> np.ndarray:
    """
    D'(y) = 0.5 * rho * V^2 * Cd(y) * c(y)
    """
    q_inf = 0.5 * rho * V_inf**2
    D_prime = q_inf * chord * ICd
    return D_prime

def compute_normal_force_distribution(L_prime: np.ndarray,
                                      aoa: float,
                                      V_inf: float = 10.0,
                                      rho: float = 1.225) -> np.ndarray:
    """
    N'(y) = L'(y) * cos(aoa) + D'(y) * sin(aoa)
    """
    aoa_rad = np.radians(aoa)
    N_prime = L_prime * np.cos(aoa_rad) + D_prime * np.sin(aoa_rad)
    return N_prime
                                      
                                      
