import sys
import os 
import math as m
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt 
from bisect import bisect_right

# ==========================================
# 1. PATH CONFIGURATION
# ==========================================
# Get the directory where this script (webbuckling.py) is located (WP5)
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

# Get the project root (one folder up from WP5)
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))

# Get the WP4 directory path
WP4_DIR = os.path.join(PROJECT_ROOT, "WP4")

# Add WP4 to sys.path
# This is crucial: it allows Integration.py to find 'XFLR', 'TorqueDist', etc.
if WP4_DIR not in sys.path:
    sys.path.insert(0, WP4_DIR)

# Add Project Root to sys.path (standard practice)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

# ==========================================
# 2. IMPORTS
# ==========================================
# We import x_grid from WP4.Integration. 
# Because we added WP4 to sys.path above, Integration.py will run without errors.
from WP4.Integration import x_grid

# ==========================================
# 3. LOADING DATA
# ==========================================
print("--- Loading Data ---")

# A. Load Spar Data (Located in WP4 folder)
# We use WP4_DIR because that is where the .npy files were generated
try:
    h_fs = np.load(os.path.join(WP4_DIR, "h_front_spar.npy"))
    h_rs = np.load(os.path.join(WP4_DIR, "h_rear_spar.npy"))
    c_upper = np.load(os.path.join(WP4_DIR, "c_upper.npy"))
    print(f"Spar Data Loaded. First values -> Front: {h_fs[0]:.4f}, Rear: {h_rs[0]:.4f}")
except FileNotFoundError:
    print(f"CRITICAL ERROR: Spar .npy files not found in {WP4_DIR}")
    print("Please run 'WP4/Integration.py' first to generate these files.")
    sys.exit(1)

# B. Load Shear & Torque Data (Located in WP4 folder)
try:
    shear  = np.load(os.path.join(WP4_DIR, "S_vals.npy"))
    torque = np.load(os.path.join(WP4_DIR, "T_torque.npy"))
    print("Shear and Torque Data Loaded.")
except FileNotFoundError:
    print(f"CRITICAL ERROR: S_vals.npy or T_torque.npy not found in {WP4_DIR}")
    sys.exit(1)

# C. Load Web Buckling CSV (Located in WP5 folder)
# We use CURRENT_DIR because this file IS inside WP5
csv_path = os.path.join(CURRENT_DIR, "ks_vs_ab.csv")
try:
    data = np.loadtxt(csv_path, delimiter=",")
    print("interpolation data loaded.")
except OSError:
    print(f"CRITICAL ERROR: Could not find '{csv_path}'. Check filename or location.")
    sys.exit(1)


# Material properties
Pois = 0.33
E = 72.4 * 10**9 # Pa
t = 0.06         # m

k_v = 1.3


# DATA PROCESSING & FUNCTIONS


ab_data = data[:, 0] * 5
ks_data = data[:, 1] * 10 + 5

_ks_interp = PchipInterpolator(ab_data, ks_data, extrapolate=False)

def average_shear_stress(V, h_s, h_r, t_f, t_r):
    return V / (h_s * t_f + h_r * t_r)

def critical_shear_stress(k_s, b):
    return m.pi ** 2 * k_s * E / (12 - 12 * Pois ** 2) * (t/b)**2

def max_shear_stress(ave_shear_stress):
    return k_v * ave_shear_stress

def ks(a_over_b):
    return _ks_interp(a_over_b)

# CALCULATIONS
# total maximum shear in a wingbox - Maximum stress due to shear stress + shear due to torsion
# for now added minus sign in front of shear so that it makes sense with torque, also A_i = h_fs * c_upper (Assumption: someone said the wingbox is a rectangle)
shear_total = - average_shear_stress(shear, h_fs, h_fs, t, t) * k_v + torque / (2 * t * h_fs * c_upper)

# minimum rib spacing - if the ribs are spaced further than this, structure will fail
rib_pos = [0] #first rib to be assumed at the root

shear_str = []
ab = []
b = []

init_pos = 0
init_ind_pos = 0
for ind_pos, pos in enumerate(x_grid):

    if ind_pos == 0:
        continue

    a_over_b = (pos - init_pos) / ((h_fs[init_ind_pos] + h_fs[ind_pos]) / 2)

    # lower bound
    if a_over_b < ab_data[0]:
        continue

    # upper bound
    if a_over_b > ab_data[-1]:
        # append at upper limit
        ab.append(ab_data[-1])
        shear_str.append(critical_shear_stress(ks(ab_data[-1]), (h_fs[init_ind_pos] + h_fs[ind_pos]) / 2))
        rib_pos.append(x_grid[ind_pos - 1])
        b.append((h_fs[init_ind_pos] + h_fs[ind_pos]) / 2)

        init_pos = x_grid[ind_pos - 1]
        init_ind_pos = ind_pos - 1
        continue 

    # inside interpolation range
    shear_crit = critical_shear_stress(ks(a_over_b), (h_fs[init_ind_pos] + h_fs[ind_pos]) / 2)

    shear_max = max(shear_total[init_ind_pos: ind_pos + 1])

    if shear_crit < shear_max:
        ab.append(a_over_b)
        shear_str.append(shear_crit)
        rib_pos.append(x_grid[ind_pos - 1])
        b.append((h_fs[init_ind_pos] + h_fs[ind_pos]) / 2)

        init_pos = x_grid[ind_pos - 1]
        init_ind_pos = ind_pos - 1


rib_pos.append(x_grid[-1]) #last rib assumed to be at the tip
#print("ab =", ab)
#print("tau =",shear_str)
#print("b =",b)
#print(rib_pos)

#ARBITRARY RIB SPACING

ribs = [0, 1,2,3,4,5,6,7,7.5,8,8.5,9,9.5 ,x_grid[-1]]

#finding the closest equivalent based on x_grid

rib_lst = []
rib_lst_index = []

for rib in ribs:
    idx = bisect_right(x_grid, rib) - 1
    idx = max(idx, 0)  # safety
    rib_lst.append(x_grid[idx])
    rib_lst_index.append(idx)

rib_lst.append(x_grid[-1])
rib_lst_index.append(len(x_grid) - 1)

#list of spar height based on rib_lst_index

h_fs_lst = [h_fs[i] for i in rib_lst_index]

x = []
y = []

for i in range(len(rib_lst) - 1):
    spacing = rib_lst[i + 1] - rib_lst[i]
    h_avg = 0.5 * (h_fs_lst[i] + h_fs_lst[i + 1])

    val = critical_shear_stress(
        ks(spacing / h_avg),
        h_avg
    )

    x.extend([rib_lst[i], rib_lst[i + 1]])
    y.extend([val, val])

# PLOTTING
plt.step(x, y, where='post')
plt.xlabel("Interval")
plt.ylabel("Function value")

plt.figure(figsize=(10, 6))

# Shear force plot
plt.subplot(3, 1, 1)
plt.plot(x_grid, shear)
plt.xlabel("x [m]")
plt.ylabel("Shear force V(x) [N]")
plt.title("Shear force distribution")
plt.grid(True)

# Torque plot
plt.subplot(3, 1, 2)
plt.plot(x_grid, torque)
plt.xlabel("x [m]")
plt.ylabel("Torque T(x) [Nm]")
plt.title("Torque distribution")
plt.grid(True)

plt.subplot(3, 1, 3)
plt.plot(x_grid, shear_total)
plt.scatter(rib_pos, [0] * len(rib_pos), color = 'r', label = 'Ribs')
plt.xlabel("x [m]")
plt.ylabel("Shear stress [Pa]")
plt.title("Maximum shear stress located in the rear spar for LC9")
plt.grid(True)



plt.tight_layout()
plt.show()