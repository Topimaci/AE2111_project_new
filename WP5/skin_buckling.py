import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# PART 1: CALCULATE APPLIED STRESS (LOAD)
# ==============================================================================

# 1. Load Data
M_vals = np.load("M_vals.npy")
x_grid = np.load("X_grid.npy")  # This is our master "z" coordinate
I_xx = np.load("I_xx.npy")
h_fs = np.load("h_front_spar.npy")
h_rs = np.load("h_rear_spar.npy")

# 2. Geometry & Neutral Axis
# Distance from top spar to neutral axis (Centroid of trapezoid)
x_c = (h_rs ** 2 + h_fs ** 2 + h_fs * h_rs) / (3 * (h_rs + h_fs))

# Distance from Neutral Axis to the skin (assuming top skin at front spar is critical)
y_norm_stress_front = h_fs - x_c

# 3. Calculate Applied Stress: sigma = M * y / I
# We use abs() because buckling cares about the magnitude of compression
stress_applied = np.abs(M_vals * y_norm_stress_front / I_xx)

# 4. Apply Cutoff (Fixing the tip singularity)
cutoff_idx = 400
if len(stress_applied) > cutoff_idx:
    cutoff_stress = stress_applied[cutoff_idx]
    stress_applied[cutoff_idx:] = cutoff_stress

# ==============================================================================
# PART 2: CALCULATE CRITICAL STRESS (STRENGTH)
# ==============================================================================

# Constants
E = 71 * 10**9     # Pa
v = 0.33           # Poisson's ratio
t = 0.005          # Skin thickness [m]

# Use the loaded x_grid as our 'z' to ensure arrays match size
z = x_grid
L = z[-1]          # Total span

# Design Parameters (Configuration 1: Fewest stringers -> Critical Case)
w_wingbox = 2.8735 * 0.3  # width of wingbox [m]

# Stringer counts
n_s_root1 = 4
n_s_tip1 = 2

n_s_root2 = 7
n_s_tip2 = 4

n_s_root3 = 9
n_s_tip3 = 5

# Calculate Spacing 'b' [m]
b_value_root = w_wingbox / (n_s_root3 + 1)
b_value_tip  = w_wingbox / (n_s_tip3 + 1)

# Create b array matching the stress array size
b = np.zeros_like(z)
k_c = np.zeros_like(z)

# Apply Distribution Logic (Split at 50% span)
midpoint = L / 2
mask_root = z < midpoint
mask_tip  = z >= midpoint

b[mask_root] = b_value_root
b[mask_tip]  = b_value_tip

k_c[:] = 4.0  # Assume simply supported

# Calculate Critical Buckling Stress
Sigma_cr_skin = (k_c * (np.pi ** 2) * E) / (12 * (1 - v ** 2)) * (t / b) ** 2

# ==============================================================================
# PART 3: MARGIN OF SAFETY & PLOTTING
# ==============================================================================

# Calculate Reserve Factor (Margin of Safety)
# RF = Strength / Load
# Note: Handle division by zero if stress is 0 (though unlikely with M_vals)
with np.errstate(divide='ignore', invalid='ignore'):
    margin_of_safety = Sigma_cr_skin / stress_applied

# Plotting
plt.figure(figsize=(12, 7))

# Plot the Margin of Safety curve
plt.plot(z, margin_of_safety, label='Reserve Factor (Buckling)', color='purple', linewidth=2)

# Add Safety Threshold Line (RF = 1.0)
plt.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='Failure Threshold (RF=1)')
plt.axhline(y=1.5, color='green', linestyle=':', label='Safety Target (RF=1.5)')

# Formatting
plt.title(f'Skin Buckling Safety Margin along Span\n(Config: Root={n_s_root3}, Tip={n_s_tip3} stringers)', fontsize=14)
plt.xlabel('Spanwise Position z [m]', fontsize=12)
plt.ylabel('Reserve Factor (Strength / Load)', fontsize=12)
plt.grid(True, which='both', linestyle='--')
plt.legend()

# Limit y-axis to make the plot readable (e.g., if root margin is huge, tip detail is lost)
# You might want to adjust this based on your results, or let it auto-scale
plt.ylim(0, max(5, np.nanmax(margin_of_safety[cutoff_idx:]))) 

plt.tight_layout()
plt.show()