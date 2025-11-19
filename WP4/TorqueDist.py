import numpy as np
import matplotlib.pyplot as plt
from XFLR import y_span, chord, Ai, Cl, ICd, Cm # Importing data from XFLR in .txt form and computing aerodynamic line load

from scipy import integrate, interpolate

# Variables

V_inf = 10.0  # Freestream velocity in m/s
rho   = 1.225 # Air density in kg/m^3
aoa_deg = 0.0   # Angle of attack in degrees


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
    

# torque distribution along the blade due to external loads or lifts etc, known as t(x)
def compute_section_moment_density(chord: np.ndarray,
                                   Cm: np.ndarray,
                                   V_inf: float,
                                   rho: float = 1.225) -> np.ndarray:
    """
    M'(y) = Cm(y) * q_inf * c(y)^2 
    """
    q_inf = 0.5 * rho * V_inf**2
    M_prime = Cm * q_inf * chord**2
    return M_prime

def build_q_d_t_functions(y_span: np.ndarray,
                      N_prime: np.ndarray,
                      M_prime: np.ndarray,
                      d0: float = 0.7):   #  m, distance from blade root to shear force center / calculation point. <---  CHANGE THIS VALUE, THIS IS A DUMMY VALUE

    mask = y_span >= 0.0
    y_half = y_span[mask]
    N_half = N_prime[mask]
    M_half = M_prime[mask]   # q(x) = N'(y)

    x = y_half
    idx = np.argsort(x)
    x_sorted = x[idx]

    q_sorted = N_half[idx]   # q(x) = N'(y)
    M_sorted = M_half[idx]   # M'(y) from Cm

    # 5. Interpolation of q(x) and d(x)
    q_func = interpolate.interp1d(
        x_sorted,
        q_sorted,
        kind="cubic",
        fill_value="extrapolate"
    )

    ### Variable d(x), from leading edge to mid wingbox

    ratio_frontspar = 0.3     # ratio of chord where front spar is located
    ratio_rearspar = 0.7      # ratio of chord where rear spar is located

    d_centroid = (ratio_frontspar + (ratio_frontspar + ratio_rearspar) / 2.0) 
    


    d_false = np.full_like(x_sorted, d0)
    d_func = interpolate.interp1d(
        x_sorted,
        d_false,
        kind="linear",
        fill_value="extrapolate"
    )

    t_func = interpolate.interp1d(
        x_sorted,
        M_sorted,
        kind="cubic",
        fill_value="extrapolate"
    )

    return x_sorted, q_func, d_func, t_func




# Torque density distribution w(x), where w(x) = q(x) * d(x)
def torque_density_distribution(x: np.ndarray,
                                q_func,
                                d_func,
                                t_func) -> np.ndarray:
    """
    w(x) = q(x) * d(x) + t(x)
    """
    w_x = q_func(x) * d_func(x) + t_func(x)
    return w_x





### Main to be completed, still test code ###
if __name__ == "__main__":

    # Check ruwe input
    print("y_span[:5] =", y_span[:5])
    print("chord[:5]  =", chord[:5])
    print("Cl[:5]     =", Cl[:5])
    print("Cm[:5]     =", Cm[:5])

    # 2. Lift and drag line loads
    L_prime = compute_lift_line_load(chord, Cl, V_inf, rho)
    D_prime = compute_drag_line_load(chord, ICd, V_inf, rho)

    print("L_prime[:5] =", L_prime[:5])
    print("D_prime[:5] =", D_prime[:5])

    # 3. Normale belasting N'(y)
    N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg)
    print("N_prime[:5] =", N_prime[:5])

    # 4. Section moment density M'(y) uit Cm
    M_prime = compute_section_moment_density(chord, Cm, V_inf, rho)
    print("M_prime[:5] =", M_prime[:5])

    # 5. q(x), d(x) en t(x) functies
    x_sorted, q_func, d_func, t_func = build_q_d_t_functions(
        y_span, N_prime, M_prime, d0=0.7
    )
    print("x_sorted[:5] =", x_sorted[:5])
    print("q(x_sorted[:5]) =", q_func(x_sorted[:5]))
    print("d(x_sorted[:5]) =", d_func(x_sorted[:5]))
    print("t(x_sorted[:5]) =", t_func(x_sorted[:5]))

    # 6. Torque-density w_T(x) op een grid
    L_span = x_sorted[-1]
    x_grid = np.linspace(0, L_span, 400)
    w_T = torque_density_distribution(x_grid, q_func, d_func, t_func=t_func)
    print("w_T[:5] =", w_T[:5])

    # 7. Integreren naar torsiediagram T(x): T(x) = -∫_x^L w_T(ξ) dξ
    x_rev = x_grid[::-1]
    w_rev = w_T[::-1]

    T_rev = integrate.cumulative_trapezoid(w_rev, x_rev, initial=0.0)
    T = -T_rev[::-1]
    print("T[:5] =", T[:5])

    # 8. Plot
    plt.figure()
    plt.plot(x_grid, T)
    plt.xlabel("Spanwise position x [m]")
    plt.ylabel("Torque T(x) [Nm]")
    plt.grid(True)
    plt.title("Torque diagram with q·d and Cm contribution")
    plt.show()
