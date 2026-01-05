import sys
import os
import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# ==========================================
# 1. PATH CONFIGURATION
# ==========================================
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
WP4_DIR = os.path.join(PROJECT_ROOT, "WP4")

for p in [WP4_DIR, PROJECT_ROOT]:
    if p not in sys.path:
        sys.path.insert(0, p)

# ==========================================
# 2. IMPORTS
# ==========================================
from WP4.Integration import x_grid

# ==========================================
# 3. LOADING DATA
# ==========================================
h_fs = np.load(os.path.join(WP4_DIR, "h_front_spar.npy"))
h_rs = np.load(os.path.join(WP4_DIR, "h_rear_spar.npy"))
c_upper = np.load(os.path.join(WP4_DIR, "c_upper.npy"))

shear  = np.load(os.path.join(WP4_DIR, "S_vals.npy"))
torque = np.load(os.path.join(WP4_DIR, "T_torque.npy"))

csv_path = os.path.join(CURRENT_DIR, "ks_vs_ab.csv")
data = np.loadtxt(csv_path, delimiter=",")

# ==========================================
# 4. MATERIAL & CONSTANTS
# ==========================================
Pois = 0.33
E = 72.4e9          # Pa
t = 0.008           # m
k_v = 1.3           # shear correction
SF = 1.0            # safety factor (increase if desired)

# ==========================================
# 5. INTERPOLATION (1 <= a/b <= 5)
# ==========================================
ab_data = data[:, 0] * 5
ks_data = data[:, 1] * 10 + 5

ks_interp = PchipInterpolator(ab_data, ks_data, extrapolate=False)

def ks(a_over_b):
    return ks_interp(a_over_b)

# ==========================================
# 6. STRESS FUNCTIONS
# ==========================================
def average_shear_stress(V, h_s, h_r, t_s, t_r):
    return V / (h_s * t_s + h_r * t_r)

def critical_shear_stress(k_s, b):
    return (m.pi**2 * k_s * E / (12 * (1 - Pois**2))) * (t / b)**2

def max_shear_stress(tau_avg):
    return k_v * tau_avg

# ==========================================
# 7. TOTAL APPLIED SHEAR STRESS
# ==========================================
tau_shear = max_shear_stress(
    average_shear_stress(shear, h_fs, h_fs, t, t)
)

tau_torsion = torque / (2 * t * h_fs * c_upper)

shear_total = np.abs(tau_shear + tau_torsion)

# ==========================================
# 8. RIB PLACEMENT ALGORITHM (CORRECT)
# ==========================================
rib_pos = [x_grid[0]]
shear_allow = []
b_vals = []
ab_vals = []

i_start = 0

while i_start < len(x_grid) - 1:

    last_safe_i = None
    last_safe_tau = None
    last_safe_b = None
    last_safe_ab = None

    for i in range(i_start + 1, len(x_grid)):

        a = x_grid[i] - x_grid[i_start]
        b_avg = 0.5 * (h_fs[i_start] + h_fs[i])
        a_over_b = a / b_avg

        # Enforce a/b limits
        if a_over_b < ab_data[0]:
            continue
        if a_over_b > ab_data[-1]:
            break

        tau_crit = critical_shear_stress(
            ks(a_over_b),
            b_avg
        )

        tau_max = np.max(shear_total[i_start:i+1])

        if tau_crit >= SF * tau_max:
            last_safe_i = i
            last_safe_tau = tau_crit
            last_safe_b = b_avg
            last_safe_ab = a_over_b
        else:
            break

    if last_safe_i is None:
        raise RuntimeError(
            f"No admissible rib spacing from x = {x_grid[i_start]:.3f} m"
        )

    rib_pos.append(x_grid[last_safe_i])
    shear_allow.append(last_safe_tau)
    b_vals.append(last_safe_b)
    ab_vals.append(last_safe_ab)

    i_start = last_safe_i

rib_pos.append(x_grid[-1])

# ==========================================
# 9. BUILD STEPWISE ALLOWABLE CURVE
# ==========================================
x_crit = []
y_crit = []

for i in range(len(shear_allow)):
    x0 = rib_pos[i]
    x1 = rib_pos[i + 1]
    x_crit.extend([x0, x1])
    y_crit.extend([shear_allow[i], shear_allow[i]])


print("Rib positions [m]:")
for r in rib_pos:
    print(f"{r:.3f}")

print("\nBay a/b ratios:")
for i, ab in enumerate(ab_vals):
    print(f"Bay {i+1}: a/b = {ab:.3f}")

# ==========================================
# 10. PLOTTING
# ==========================================
plt.figure(figsize=(10, 8))

plt.subplot(3, 1, 1)
plt.plot(x_grid, shear)
plt.ylabel("Shear force [N]")
plt.grid(True)
plt.title("Shear force distribution")

plt.subplot(3, 1, 2)
plt.plot(x_grid, torque)
plt.ylabel("Torque [Nm]")
plt.grid(True)
plt.title("Torque distribution")

plt.subplot(3, 1, 3)
plt.plot(x_grid, shear_total, label="Applied shear stress")
plt.step(x_crit, y_crit, where="post", color="orange",
         label="Allowable shear stress")
plt.scatter(rib_pos, [0]*len(rib_pos), color="red", label="Ribs")
plt.xlabel("x [m]")
plt.ylabel("Shear stress [Pa]")
plt.legend()
plt.grid(True)
plt.title("Shear buckling–driven rib spacing")

plt.tight_layout()
plt.show()