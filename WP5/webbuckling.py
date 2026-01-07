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
from WP4.Stiffness import t_spar

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
t = t_spar           # m 
k_v = 1.3           # shear correction
SF = 1.0            # safety factor

# ==========================================
# 5. INTERPOLATION (1 <= a/b <= 5, roughly not exactly 1 and 5)
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
# 8. RIB PLACEMENT ALGORITHM for finding the maximum rib spacing before failure
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

#====
#10. Arbitrary rib spacing - FINAL RIB SPACING
#====

#If you get a/b is larger error, just adjust the rib spacing accordingly
rib_position = [0, 1.685, 3.174, 4.506, 5.681, 6.739 , 7.68, 8.522, 9.266, x_grid[-1]]

rib_pos_new = []
rib_pos_new_ind = []

# Finding closest LOWER x_grid positions
last_ind = -1

for pos in rib_position:
    ind = int(np.argmin(np.abs(x_grid - pos)))

    # enforce strict monotonicity
    if ind <= last_ind:
        continue

    rib_pos_new_ind.append(ind)
    rib_pos_new.append(x_grid[ind])
    last_ind = ind

x_vals = []
y_vals = []

for i in range(len(rib_pos_new) - 1):
    rib1 = rib_pos_new[i]
    rib2 = rib_pos_new[i + 1]

    height1 = h_fs[rib_pos_new_ind[i]]
    height2 = h_fs[rib_pos_new_ind[i + 1]]

    avg_height = 0.5 * (height1 + height2)
    spacing = rib2 - rib1

    ab_val = spacing / avg_height
    k_s = ks(ab_val)
    #Making sure a/b is in range
    if ab_val > ab_data[-1]:
        print("The a/b ratio is larger than 5 between ribs at ", rib1, " and ", rib2)
        break
    if ab_val < ab_data[0]:
        print("The a/b ratio is smaller than 1 between ribs at ", rib1, " and ", rib2)
        break

    cr_shear = critical_shear_stress(k_s, avg_height)

    x_vals.extend([rib1, rib2])
    y_vals.extend([cr_shear, cr_shear])

avg_bay_stress = []

for i in range(len(rib_pos_new_ind) - 1):
    i1 = rib_pos_new_ind[i]
    i2 = rib_pos_new_ind[i + 1]

    # Average applied shear stress in bay
    tau_avg_bay = np.max(shear_total[i1:i2+1])

    avg_bay_stress.extend([tau_avg_bay, tau_avg_bay])
# ==========================================
# 11. PLOTTING
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
#Comment / uncomment if you need the maximum rib spacing graph
plt.step(x_crit, y_crit, where="post", color="orange", label="Allowable shear stress")
plt.scatter(rib_pos, [0]*len(rib_pos), color="red", label="Ribs") 
plt.step(x_vals, y_vals, where="post", label = "Shear - chosen ribs", color="green")
plt.scatter(rib_pos_new, [0]*len(rib_pos_new), color="blue", label="Ribs final")
plt.xlabel("x [m]")
plt.ylabel("Shear stress [Pa]")
plt.legend()
plt.grid(True)
plt.title("Shear buckling–driven rib spacing")

plt.tight_layout()
plt.show()


#FOR ORIOL
safety_margin_x_values = x_vals
margin_of_safety = np.array(y_vals) / np.array(avg_bay_stress)

print(margin_of_safety)
plt.step(x_vals, margin_of_safety)
plt.show()