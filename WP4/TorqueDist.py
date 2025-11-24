import numpy as np
import matplotlib.pyplot as plt
from XFLR import y_span, chord, Ai, Cl, ICd, Cm # Importing data from XFLR in .txt form and computing aerodynamic line load,
import math as m
import math as m
from matplotlib.widgets import RadioButtons
from scipy import integrate, interpolate

# Variables

V_inf = 52.74  # Freestream velocity in m/s
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
                      d0: float = 0.7):   #  m, distance from load to flexural axis / calculation point. <---  CHANGE THIS VALUE, THIS IS A DUMMY VALUE

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

def distance_dx_calc(chord, Cl, Cm):
    # distance = 0.45*c − 0.25*c = 0.20*c
    # 0.45 comes from middle of spars, 20% and 70% needs to be checked (0.45c from the LE)
    # 0.25 comes from assumption that lift acts as a point force on the 0.25 c from the LE
    #d_extra = Cm/Cl
    dx = 0.45*chord - 0.25*chord 
    sweep_deg = 10.43 #from WP3 sweep at quarter chord
    sweep_rad = m.radians(sweep_deg)
    dreal = dx * m.cos(sweep_rad)
    return dreal





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

# Extra forces on the wing, e.g. fuel weight, engine weight, etc. can be added here as additional torque densities if needed.
def add_point_forces_and_torques(x_grid: np.ndarray,
                                 T_dist: np.ndarray,
                                 point_forces=None,
                                 point_torques=None):
    """
    Add point forces and torques to the torque distribution.
    point_forces: list of tuples (position, magnitude)
    point_torques: list of tuples (position, magnitude)
    """ 
    T_total = T_dist.copy()

    if point_forces is not None:
        for pf in point_forces:
            xP = pf['x']
            P = pf['P']
            d = pf['d']
            T_total += P * d * (x_grid <= xP)

    if point_torques is not None:
        for pt in point_torques:
            xT = pt['x']
            T_mag = pt['T']
            T_total += T_mag * (x_grid <= xT)

    return T_total



if __name__ == "__main__":

    # === your existing computations ===
    print("y_span[:5] =", y_span[:5])
    print("chord[:5]  =", chord[:5])
    print("Cl[:5]     =", Cl[:5])
    print("Cm[:5]     =", Cm[:5])

    L_prime = compute_lift_line_load(chord, Cl, V_inf, rho)
    D_prime = compute_drag_line_load(chord, ICd, V_inf, rho)
    N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg)
    M_prime = compute_section_moment_density(chord, Cm, V_inf, rho)

    x_sorted, q_func, d_func, t_func = build_q_d_t_functions(
        y_span, N_prime, M_prime, d0=0.7
    )

    L_span = x_sorted[-1]
    x_grid = np.linspace(0, L_span, 400)
    w_T = torque_density_distribution(x_grid, q_func, d_func, t_func=t_func)

    x_rev = x_grid[::-1]
    w_rev = w_T[::-1]
    T_rev = integrate.cumulative_trapezoid(w_rev, x_rev, initial=0.0)
    T_dist = -T_rev[::-1]

    point_forces = [
        {'x': 1.84, 'P': 126.8*9.81, 'd': 0.473},  # Landing gear

    ]

    point_torques = [
        #{'x': 1.84, 'T': -point_forces[0]['P']*point_forces[0]['d']},  # torque due to landing gear weight
        {'x': 1.84, 'T': 0.5 * rho * V_inf**2*0.04905*(0.56/2)**2*m.pi*0.785}  # torque due to landing gear drag
    
    ]
    print(point_torques)


    T_total = add_point_forces_and_torques(
        x_grid, T_dist,
        point_forces=point_forces,
        point_torques=point_torques
    )

    # === INTERACTIVE UI WITH RADIO BUTTONS ===

    # 1) Make figure + axis, leave room for the radio buttons on the left
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.3)

    # 2) Initial plot: show the torque diagram
    def plot_torque():
        ax.clear()
        ax.plot(x_grid, T_dist, label="Distributed loads only")
        ax.plot(x_grid, T_total, label="With point forces/torques")
        ax.set_xlabel("Spanwise position x [m]")
        ax.set_ylabel("Torque T(x) [Nm]")
        ax.set_title("Torque diagram")
        ax.grid(True)
        ax.legend()

    def plot_line_loads():
        ax.clear()
        ax.plot(y_span, L_prime, label="Lift line load L'(y)")
        ax.plot(y_span, D_prime, label="Drag line load D'(y)")
        ax.set_xlabel("Spanwise position y [m]")
        ax.set_ylabel("Line load [N/m]")
        ax.set_title("Aerodynamic line loads")
        ax.grid(True)
        ax.legend()

    # Start with the torque view
    plot_torque()

    # 3) Radio buttons on the left
    ax_radio = plt.axes([0.05, 0.4, 0.2, 0.15])  # [left, bottom, width, height]
    labels = ["Torque", "Line loads"]
    radio = RadioButtons(ax_radio, labels)

    # 4) Callback: what happens when you click a radio button
    def change_plot(label):
        if label == "Torque":
            plot_torque()
        elif label == "Line loads":
            plot_line_loads()
        fig.canvas.draw_idle()

    radio.on_clicked(change_plot)

    # 5) Show the UI
    plt.show()

