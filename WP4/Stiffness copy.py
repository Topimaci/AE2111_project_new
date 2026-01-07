import math as m
import numpy as np
from scipy.interpolate import interp1d
from Integration import x_grid
from Moment_Diagram import M_vals

# --- Constants ---
b = 19.585
root_chord = 2.874
tip_chord = 1.043

# --- Airfoil Data ---
Airfoil_coordinates = []
with open("WP4/NACA64714 a=0.0.dat", "r") as file:
    for line in file:
        line = line.strip()
        if not line or not line[0].isdigit():
            continue
        x, y = map(float, line.split())
        Airfoil_coordinates.append([x, y])

# --- Geometry Helper Functions ---
def spar_position(coords_list, fraction):
    c1 = c2 = c3 = c4 = [0, 0]
    prev = [0, 0]
    one_found = True
    for coords in coords_list:
        if fraction > coords[0] and one_found:
            c2, c1, one_found = coords, prev, False
        if fraction < coords[0] and not one_found:
            c4, c3 = coords, prev
            break
        prev = coords
    return c1, c2, c3, c4

def top_stringer_y_coord(s1c1, s1c2, s2c1, s2c2, frac1):
    y_top1, y_top2 = max(s1c1[1], s1c2[1]), max(s2c1[1], s2c2[1])
    top_y = min(y_top1, y_top2)
    y_bot1, y_bot2 = min(s1c1[1], s1c2[1]), min(s2c1[1], s2c2[1])
    bot_y = max(y_bot1, y_bot2)
    if top_y - bot_y <= 0:
        return 0.5 * (y_top1 + y_top2)
    a = (s1c2[1] - s1c1[1]) / (s1c2[0] - s1c1[0])
    return a * frac1 + (s1c1[1] - a * s1c1[0])

def bot_stringer_y_coords(s1c3, s1c4, s2c3, s2c4, frac1, frac2):
    a1 = (s1c4[1] - s1c3[1]) / (s1c4[0] - s1c3[0])
    y1 = a1 * frac1 + (s1c3[1] - a1 * s1c3[0])
    a2 = (s2c4[1] - s2c3[1]) / (s2c4[0] - s2c3[0])
    y2 = a2 * frac2 + (s2c3[1] - a2 * s2c3[0])
    return y1, y2

def spar_length(frac, y_pos, r_c, t_c, span, f1, f2):
    s1c1, s1c2, s1c3, s1c4 = spar_position(Airfoil_coordinates, f1)
    s2c1, s2c2, s2c3, s2c4 = spar_position(Airfoil_coordinates, f2)
    y_top = top_stringer_y_coord(s1c1, s1c2, s2c1, s2c2, f1)
    yb1, yb2 = bot_stringer_y_coords(s1c3, s1c4, s2c3, s2c4, f1, f2)
    chord_y = r_c + ((t_c - r_c) / (span / 2)) * y_pos
    a_box = (yb2 - yb1) / (f2 - f1)
    y_bot_at_spar = a_box * frac + (yb1 - a_box * f1)
    return (y_top - y_bot_at_spar) * chord_y

def stiffness_distribution(y_pos, h_fs, h_rs, c_up, c_low, t_sk, t_sp, A_str, interp_top, interp_bot):
    x_c = (h_rs**2 + h_fs**2 + h_fs*h_rs) / (3*(h_rs + h_fs))
    I_fs = (1/12 * h_fs**3 * t_sp) + (h_fs * t_sp * (x_c - h_fs/2)**2)
    I_rs = (1/12 * h_rs**3 * t_sp) + (h_rs * t_sp * (x_c - h_rs/2)**2)
    I_top = t_sk * c_up * (t_sk/2 - x_c)**2
    I_bot = t_sk * c_low * (((h_fs - x_c) + (h_rs - x_c))/2)**2
    num_t, num_b = max(0, interp_top(y_pos)), max(0, interp_bot(y_pos))
    I_str_t = (A_str * (t_sk - x_c)**2) * num_t
    I_str_b = (A_str * (((h_fs - x_c) + (h_rs - x_c))/2)**2) * num_b
    return I_fs + I_rs + I_top + I_bot + I_str_t + I_str_b

# --- Design Definitions ---
designs = {
    "D1": {"str_t": [4,4,2,2], "str_b":[4,4,2,2], "t_sk":0.002, "t_sp":0.005, "A_str":0.0001},
    "D2": {"str_t": [7,7,4,4], "str_b":[7,7,4,4], "t_sk":0.003, "t_sp":0.008, "A_str":0.00025},
    "D3": {"str_t": [9,9,5,5], "str_b":[9,9,5,5], "t_sk":0.003, "t_sp":0.005, "A_str":0.0002},
    "D4": {"str_t": [8, 7, 5, 4], "str_b":[5, 5, 4, 4], "t_sk":0.005, "t_sp":0.010, "A_str":0.0004},
    "D5": {"str_t": [9, 8, 6, 4], "str_b":[7, 7, 4, 3], "t_sk":0.004, "t_sp":0.008, "A_str":0.0003},
}

y_breaks = np.array([0, 3, 4.89, 7])

for name, d in designs.items():
    itop = interp1d(y_breaks, d["str_t"], kind='linear', fill_value="extrapolate")
    ibot = interp1d(y_breaks, d["str_b"], kind='linear', fill_value="extrapolate")
    
    I_vals = []
    h_fs_arr = []
    h_rs_arr = []
    
    for yi in x_grid:
        h_fs = spar_length(0.3, yi, root_chord, tip_chord, b, 0.3, 0.6)
        h_rs = spar_length(0.6, yi, root_chord, tip_chord, b, 0.3, 0.6)
        chord_y = root_chord + ((tip_chord - root_chord)/(b/2)) * yi
        c_up = abs(0.6 - 0.3) * chord_y
        c_low = c_up  # simplified
        I = stiffness_distribution(yi, h_fs, h_rs, c_up, c_low, d["t_sk"], d["t_sp"], d["A_str"], itop, ibot)
        I_vals.append(I)
        h_fs_arr.append(h_fs)
        h_rs_arr.append(h_rs)
    
    # Save everything needed for MoS calculation
    np.save(f"I_xx_{name}.npy", np.array(I_vals))
    np.save(f"h_front_spar_{name}.npy", np.array(h_fs_arr))
    np.save(f"h_rear_spar_{name}.npy", np.array(h_rs_arr))
    print(f"Saved I_xx, h_front_spar, h_rear_spar for {name}")
