import sys
import os 
import math as m
import numpy as np
from scipy.interpolate import PchipInterpolator
import matplotlib.pyplot as plt 

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
    print("Web buckling CSV loaded.")
except OSError:
    print(f"CRITICAL ERROR: Could not find '{csv_path}'. Check filename or location.")
    sys.exit(1)

# ==========================================
# 4. DATA PROCESSING & FUNCTIONS
# ==========================================

ab_data = data[:, 0] * 5
ks_data = data[:, 1] * 10 + 5

_ks_interp = PchipInterpolator(ab_data, ks_data, extrapolate=False)

# Material properties
Pois = 0.33
E = 72.4 * 10**9 # Pa
t = 0.06         # m

k_v = 1.3

def average_shear_stress(V, h_s, h_r, t_f, t_r):
    return V / (h_s * t_f + h_r * t_r)

def critical_shear_stress(k_s, b):
    # Note: I added **2 to (t/b) as that is the standard buckling formula.
    # If your project explicitly requires linear (t/b), remove the **2.
    return m.pi ** 2 * k_s * E / (12 - 12 * Pois ** 2) * (t/b)**2

def max_shear_stress(ave_shear_stress):
    return k_v * ave_shear_stress

def ks(a_over_b):
    return _ks_interp(a_over_b)

# ==========================================
# 5. PLOTTING
# ==========================================
plt.figure(figsize=(10, 6))

# Shear force plot
plt.subplot(2, 1, 1)
plt.plot(x_grid, shear)
plt.xlabel("x [m]")
plt.ylabel("Shear force V(x) [N]")
plt.title("Shear force distribution")
plt.grid(True)

# Torque plot
plt.subplot(2, 1, 2)
plt.plot(x_grid, torque)
plt.xlabel("x [m]")
plt.ylabel("Torque T(x) [Nm]")
plt.title("Torque distribution")
plt.grid(True)

plt.tight_layout()
plt.show()