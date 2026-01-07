import numpy as np
import matplotlib.pyplot as plt

# --- SECTION 1: LOAD DATA & CALCULATE APPLIED STRESS (Your existing code) ---
M_vals = np.load("M_vals.npy")
x_grid = np.load("X_grid.npy")
I_xx = np.load("I_xx.npy")
h_fs = np.load("h_front_spar.npy")
h_rs = np.load("h_rear_spar.npy")

# Distance from top spar to neutral axis
x_c = (h_rs ** 2 + h_fs ** 2 + h_fs * h_rs) / (3 * (h_rs + h_fs))

# Distance from neutral axis to the skin (Critical fiber distance)
y_norm_stress_front = h_fs - x_c

# Calculate Applied Stress
stress = M_vals * y_norm_stress_front / I_xx

# Apply the cutoff (Tip correction)
cutoff_value = 400
# Ensure we don't go out of bounds if array is smaller than 500
if len(stress) > cutoff_value:
    cutoff_stress = stress[cutoff_value]
    stress[cutoff_value:] = cutoff_stress

# --- SECTION 2: DEFINE GEOMETRY & CONSTANTS (From your first code) ---
# We re-define these here to ensure they match the size of 'x_grid'
E = 71 * 10**9   # Pa
v = 0.33         # Poisson's ratio
t = 0.005        # Skin thickness [m]

# Re-create the 'b' array using 'x_grid' instead of a new linspace
# This ensures 'b' and 'stress' arrays are exactly the same length
b = np.zeros_like(x_grid)

# Your stringer spacing inputs
w_wingbox = 2.8735 * 0.3 
n_s_root3 = 9
n_s_tip3  = 5

# Corrected Spacing Formulas (Meters)
b_value_root = w_wingbox / (n_s_root3 + 1)
b_value_tip  = w_wingbox / (n_s_tip3 + 1)

# Apply spacing logic based on x_grid
L = x_grid[-1] # Assume total length is the last point in x_grid
midpoint = L / 2

b[x_grid < midpoint]  = b_value_root
b[x_grid >= midpoint] = b_value_tip

# --- SECTION 3: SOLVE FOR REQUIRED k_c ---

# Formula derived: k_c = (Sigma_applied * 12 * (1-v^2)) / (pi^2 * E * (t/b)^2)
# We use np.abs(stress) to ensure we are looking at magnitude (compressive load)
term1 = 12 * (1 - v**2)
term2 = np.pi**2 * E * (t / b)**2
k_c_needed = (np.abs(stress) * term1) / term2

# --- SECTION 4: OUTPUT & PLOTTING ---

print("Max k_c required:", np.max(k_c_needed))
print("Mean k_c required:", np.mean(k_c_needed))

plt.figure(figsize=(10, 6))
plt.plot(x_grid, k_c_needed, label='Required $k_c$', color='purple', linewidth=2)

# Add a reference line for k_c = 4 (Standard simply supported skin)
plt.axhline(y=4, color='green', linestyle='--', label='Standard limit ($k_c=4$)')

plt.xlabel('Spanwise Location [m]')
plt.ylabel('Buckling Coefficient $k_c$')
plt.title('Required Skin Buckling Coefficient ($k_c$) to prevent Failure')
plt.grid(True, which="both", linestyle='--')
plt.legend()
plt.show()