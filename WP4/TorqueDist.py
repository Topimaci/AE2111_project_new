import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")    #so graph can be interactive in pycharm
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons, Slider
from XFLR import (
    y_span0, chord0, Ai0, Cl0, ICd0, Cm0,
    y_span10, chord10, Ai10, Cl10, ICd10, Cm10
)
import math as m
from matplotlib.widgets import RadioButtons
from scipy import integrate, interpolate
from shear_centre_location import shear_center_non_dim
from data_for_weight_loads_torsion import combined_loads_weights_wing_fuel


# Variables

V_inf = 200  # Freestream velocity in m/s
rho   = 1.225 # Air density in kg/m^3
   # Angle of attack in degrees

 # Angle of attack in degrees
n = 2.68
S_wing = 38.379
W_situation = 140000     #### FORCE NOT MASS
#OEM: 7881 = 77312 N
###MTOM: 14266 = 140000 N
###Payload design: 750kg =7357 N


### code for load factor and determining critical alpha
def critical_alpha(rho, v_situation, S_wing, W_situation, n, landing = False, takeoff = False):
    CL = 2*n*W_situation/(rho*v_situation**2*S_wing)
    print(CL)
    if landing == True:
        CL -= 1.15
    if takeoff == True:
        CL -= 1.15
    ### from simulation
    aoa_critical = (CL-0.327220)*10/(1.218686-0.327220)     

    return aoa_critical    


aoa_deg = critical_alpha(rho, V_inf, S_wing, W_situation, n)

print(aoa_deg)



def compute_lift_line_load(chord: np.ndarray,
                           Cl0: np.ndarray,
                           Cl10: np.ndarray,
                            aoa_deg: float,
                           V_inf: float,
                           rho: float = 1.225) -> np.ndarray:
    """
    L'(y) = 0.5 * rho * V^2 * Cl(y) * c(y)
    """

    Cl = (Cl10 - Cl0) / 10 * aoa_deg + Cl0
    

    q_inf = 0.5 * rho * V_inf**2
    L_prime = q_inf * chord * Cl
    return L_prime

def compute_drag_line_load(chord: np.ndarray, ICd0: np.ndarray, ICd10: np.ndarray, aoa_deg: float, V_inf: float, rho: float = 1.225) -> np.ndarray:
    """
    D'(y) = 0.5 * rho * V^2 * Cd(y) * c(y)
    """
    ICd = (ICd10 - ICd0) / 10 * aoa_deg + ICd0

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
                                   Cm0: np.ndarray,
                                   Cm10: np.ndarray,
                                   aoa_deg: float,
                                   V_inf: float,
                                   rho: float = 1.225) -> np.ndarray:
    """
    M'(y) = Cm(y) * q_inf * c(y)^2 
    """

    Cm = (Cm10 - Cm0) / 10 * aoa_deg + Cm0

    q_inf = 0.5 * rho * V_inf**2
    M_prime = Cm * q_inf * chord**2
    return M_prime

def build_q_d_t_functions(y_span: np.ndarray,
                          chord: np.ndarray,
                          N_prime: np.ndarray,
                          M_prime: np.ndarray,
                          sweep_deg: float = 10.43,
                          ratio_frontspar: float = 0.3,
                          ratio_rearspar: float = 0.7,
                          x_force_ratio: float = 0.25):

    mask = y_span >= 0.0
    y_half = y_span[mask]
    N_half = N_prime[mask]
    M_half = M_prime[mask]
    c_half = chord[mask]

    x = y_half
    idx = np.argsort(x)
    x_sorted = x[idx]

    q_sorted = N_half[idx]
    M_sorted = M_half[idx]
    c_sorted = c_half[idx]

    q_func = interpolate.interp1d(x_sorted, q_sorted, kind="cubic",  fill_value="extrapolate")
    t_func = interpolate.interp1d(x_sorted, M_sorted, kind="cubic",  fill_value="extrapolate")

    
    dreal_sorted = distance_dx_calc(
        c_sorted,
        ratio_frontspar=ratio_frontspar,
        ratio_rearspar=ratio_rearspar,
        x_force_ratio=x_force_ratio,
        sweep_deg=sweep_deg
    )

    d_func = interpolate.interp1d(x_sorted, dreal_sorted, kind="linear", fill_value="extrapolate")

    return x_sorted, q_func, d_func, t_func



def distance_dx_calc(chord: np.ndarray,
                     ratio_frontspar: float = 0.3,
                     ratio_rearspar: float = 0.7,
                     x_force_ratio: float = 0.25,
                     sweep_deg: float = 10.43) -> np.ndarray:
    """
    dreal(x) = (x_wingbox - x_force) * cos(sweep)
    x_wingbox = average position of front and rear spar
    x_force   = position of aerodynamic force (25% chord)
    """
    chord = np.asarray(chord)

    x_wb   = shear_center_non_dim() * chord
    x_force = x_force_ratio * chord

    sweep_rad = np.deg2rad(sweep_deg)
    dreal = (x_wb - x_force) * np.cos(sweep_rad)
    return dreal

def distance_dx_calc_wing_load_distribution(chord: np.ndarray,
                     ratio_frontspar: float = 0.3,
                     ratio_rearspar: float = 0.7,
                     x_force_ratio: float = 0.45,
                     sweep_deg: float = 8.36) -> np.ndarray:
    """
    dreal(x) = (SC - x_force) * cos(sweep)
    x_force   = centroid of wing-box (45% chord)
    """
    chord = np.asarray(chord)

    x_wb   = shear_center_non_dim() * chord
    x_force = x_force_ratio * chord

    sweep_rad = np.deg2rad(sweep_deg)
    dreal = (x_wb - x_force) * np.cos(sweep_rad)
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

    safety_zone = (x_grid >= 0.0) & (x_grid <= 0.78)
    T_total[safety_zone] *= 2.8


    return T_total


def compute_case(y_span, chord, Cl0, Cl10, aoa_deg, ICd0, ICd10, Cm0, Cm10, V_inf, rho):
    # 1. Lift & drag
    L_prime = compute_lift_line_load(chord, Cl0, Cl10, aoa_deg, V_inf, rho)
    D_prime = compute_drag_line_load(chord, ICd0, ICd10, aoa_deg, V_inf, rho)

    L_total = total_from_line_load(y_span, L_prime)
    D_total = total_from_line_load(y_span, D_prime)
    #BLYET print(f"AoA={aoa_deg_case:>4.1f}°  Lift={L_total:,.1f} N   Drag={D_total:,.1f} N")

    # 2. Normal force
    N_prime = compute_normal_force_distribution(L_prime, D_prime, aoa_deg)

    # 3. Section moment density
    M_prime = compute_section_moment_density(chord, Cm0, Cm10, aoa_deg, V_inf, rho)

    # 4. q(x), d(x), t(x)
    x_sorted, q_func, d_func, t_func = build_q_d_t_functions(
        y_span=y_span, chord=chord, N_prime=N_prime, M_prime=M_prime,
        sweep_deg=10.43, ratio_frontspar=0.3, ratio_rearspar=0.7, x_force_ratio=0.25
    )

    # 5. Torque density & torsion diagram:
    L_span = x_sorted[-1]
    x_grid = np.linspace(0, L_span, 500)
    w_T = torque_density_distribution(x_grid, q_func, d_func, t_func=t_func)
    x_rev = x_grid[::-1]
    w_rev = w_T[::-1]
    T_rev = integrate.cumulative_trapezoid(w_rev, x_rev, initial=0.0)
    T_dist = -T_rev[::-1]



    # 6. Torque from weights
    Nw = combined_loads_weights_wing_fuel.size
    y_w = np.linspace(0.0, L_span, Nw)          # L_span = half-span in meters

    w_on_grid = np.interp(x_grid, y_w, combined_loads_weights_wing_fuel)  # (500,)

    chord_half = chord[y_span >= 0.0]
    y_half = y_span[y_span >= 0.0]
    chord_on_grid = np.interp(x_grid, y_half, chord_half)
    d_weight_on_grid = distance_dx_calc_wing_load_distribution(
        chord=chord_on_grid,
        x_force_ratio=0.45,
        sweep_deg=8.36
    )
    T_dist += w_on_grid * d_weight_on_grid

    

    # point loads...
    point_forces = [{'x': 1.84, 'P': 126.8*9.81, 'd': 0.602}]
    point_torques = [{'x': 1.84, 'T': 0.5 * rho * V_inf**2*0.04905*0.233*0.689624}]
    

    T_total = add_point_forces_and_torques(x_grid, T_dist, point_forces, point_torques)

    return {
        "y_span": y_span,
        "L_prime": L_prime,
        "D_prime": D_prime,
        "x_grid": x_grid,
        "T_dist": T_dist,
        "T_total": T_total,
        "L_total": L_total,
        "D_total": D_total,
    }

    



def total_from_line_load(y, fprime):
    y = np.asarray(y); fprime = np.asarray(fprime)
    idx = np.argsort(y)
    y_s = y[idx]; f_s = fprime[idx]

    # Remove duplicates
    y_u, unique_idx = np.unique(y_s, return_index=True)
    f_u = f_s[unique_idx]

    return np.trapezoid(f_u, y_u)  # N


if __name__ == "__main__":

    current_plot_type = "Torque"

    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.35, bottom=0.15)

    # --- Slider: AoA --- <---- How beautifyl is this slider???
    ax_aoa = plt.axes((0.35, 0.05, 0.6, 0.03))
    aoa_slider = Slider(ax_aoa, "AoA [deg]", 0.0, 10.0, valinit=aoa_deg, valstep=0.05)

    def update_plot(_=None):
        ax.clear()
        aoa_now = aoa_slider.val

        res = compute_case(
            y_span0, chord0,
            Cl0, Cl10,
            aoa_now,
            ICd0, ICd10,
            Cm0, Cm10,
            V_inf, rho
        )

        if current_plot_type == "Torque":
            ax.plot(res["x_grid"], res["T_dist"], label="Distributed loads only")
            ax.plot(res["x_grid"], res["T_total"], label="With point forces/torques")
            ax.set_xlabel("Spanwise position x [m]")
            ax.set_ylabel("Torque T(x) [Nm]")
            ax.set_title(f"Torque diagram (AoA = {aoa_now:.2f}°)")
        else:
            ax.plot(res["y_span"], res["L_prime"], label="Lift line load L'(y)")
            ax.plot(res["y_span"], res["D_prime"], label="Drag line load D'(y)")
            ax.set_xlabel("Spanwise position y [m]")
            ax.set_ylabel("Line load [N/m]")
            ax.set_title(f"Aerodynamic line loads (AoA = {aoa_now:.2f}°)")

        ax.grid(True)
        ax.legend()
        fig.canvas.draw_idle()

    # Slider updates
    aoa_slider.on_changed(update_plot)

    # --- Radio buttons: plot type, switch between load and torque ---
    ax_radio_plot = plt.axes((0.05, 0.65, 0.25, 0.25))
    radio_plot = RadioButtons(ax_radio_plot, ["Torque", "Line loads"])

    def on_plot_change(label):
        nonlocal_plot_type = label  # avoid "global" mess? not possible in this scope

    def on_plot_change(label):
        global current_plot_type
        current_plot_type = label
        update_plot()

    radio_plot.on_clicked(on_plot_change)

    # initial draw
    update_plot()

    plt.show()








