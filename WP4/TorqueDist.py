import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scipy import integrate, interpolate

# Variables

V_inf = 10.0  # Freestream velocity in m/s
rho   = 1.225 # Air density in kg/m^3
aoa_deg = 0.0   # Angle of attack in degrees

# Importing data from XFLR in .txt form and computing aerodynamic line load

filename = 'angleofattack0.txt'


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
                                      D_prime: np.ndarray,
                                      aoa_deg: float) -> np.ndarray:
    """
    N'(y) = L'(y) * cos(aoa) + D'(y) * sin(aoa)
    """
    aoa_rad = np.radians(aoa_deg)
    N_prime = L_prime * np.cos(aoa_rad) + D_prime * np.sin(aoa_rad)
    return N_prime # Normal force distribution N'(y) = q(x)

def functions_q_and_d(y_span: np.ndarray,
                      N_prime: np.ndarray,
                      d0: float = 0.7):   #  m, distance from blade root to shear force center / calculation point. <---  CHANGE THIS VALUE, THIS IS A DUMMY VALUE

    x = np.abs(y_span)
    idx = np.argsort(x)
    x_sorted = x[idx]
    q_sorted = N_prime[idx]   # q(x) = N'(y)

    # 5. Interpolation of q(x) and d(x)
    q_func = interpolate.interp1d(
        x_sorted,
        q_sorted,
        kind="cubic",
        fill_value="extrapolate"
    )

    d_false = np.full_like(x_sorted, d0)
    d_func = interpolate.interp1d(
        x_sorted,
        d_false,
        kind="linear",
        fill_value="extrapolate"
    )

    return q_func, d_func, x_sorted


# Torque density distribution w(x), where w(x) = q(x) * d(x)
def torque_density_distribution(x: np.ndarray,
                                q_func,
                                d_func) -> np.ndarray:
    """
    w(x) = q(x) * d(x)
    """
    w_x = q_func(x) * d_func(x)
    return w_x
    

# torque distribution along the blade due to external loads or lifts etc, known as t(x)
def distributed_torque()



# Torque to be integrated along the blade: w_T = w(x) + t(x)

if __name__ == "__main__":
    # 1. Data 
    y_span, chord, Ai, Cl, ICd, Cm = load_xflr_data(filename)

    # 2. Lift and drag line loads
    L_prime = compute_lift_line_load(chord, Cl, V_inf, rho)
    D_prime = compute_drag_line_load(chord, ICd, V_inf, rho)

    # 3. N'(x) = q(x)
    N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg)

    # 4. q(x) and d(x) functions
    q_func, d_func, x_sorted = functions_q_and_d(y_span, N_prime, d0=0.7)

    # 5. Torque density distribution w(x)

