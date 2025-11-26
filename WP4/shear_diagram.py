###Import packages
import numpy as np


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



