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
# PART 2 & 3: CALCULATE CRITICAL STRESS & PLOT FOR ALL DESIGNS
# ==============================================================================

# Constants
E = 71 * 10**9     # Pa
v = 0.33           # Poisson's ratio
z = x_grid
L = z[-1]          # Total span
w_wingbox = 2.8735 * 0.3  # width of wingbox [m]

# --- DEFINE THE 3 DESIGNS ---
# Format: Label, thickness (t), n_root, n_tip, color
designs = [
    {
        "label": "Design 1", 
        "t": 0.002, 
        "n_root": 4, 
        "n_tip": 2, 
        "color": "blue"
    },
    {
        "label": "Design 2", 
        "t": 0.006, 
        "n_root": 7, 
        "n_tip": 4, 
        "color": "orange"
    },
    {
        "label": "Design 3", 
        "t": 0.003, 
        "n_root": 9, 
        "n_tip": 5, 
        "color": "green"
    }
]

# Initialize Plot
plt.figure(figsize=(12, 7))

# Loop through each design to calculate and plot
for d in designs:
    t_current = d["t"]
    n_root = d["n_root"]
    n_tip = d["n_tip"]
    
    # Calculate Spacing 'b' [m] for this design
    b_value_root = w_wingbox / (n_root + 1)
    b_value_tip  = w_wingbox / (n_tip + 1)

    # Create b array matching the stress array size
    b = np.zeros_like(z)
    
    # Apply Distribution Logic (Split at 50% span)
    midpoint = L / 2
    mask_root = z < midpoint
    mask_tip  = z >= midpoint

    b[mask_root] = b_value_root
    b[mask_tip]  = b_value_tip

    # Buckling Coefficient
    k_c = 4.0  # Assume simply supported (vectorized or scalar works)

    # Calculate Critical Buckling Stress for this design
    Sigma_cr_skin = (k_c * (np.pi ** 2) * E) / (12 * (1 - v ** 2)) * (t_current / b) ** 2

    # Calculate Reserve Factor (Margin of Safety)
    with np.errstate(divide='ignore', invalid='ignore'):
        margin_of_safety = Sigma_cr_skin / stress_applied

    # Plot this design's curve
    plt.plot(z, margin_of_safety, label=f'{d["label"]} (t={t_current*1000}mm)', color=d["color"], linewidth=2)


# --- FORMATTING THE GRAPH ---

# Add Safety Threshold Lines
plt.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='Failure Threshold (RF=1)')
plt.axhline(y=1.5, color='black', linestyle=':', label='Safety Target (RF=1.5)')

plt.title('Skin Buckling Safety Margin Comparison', fontsize=14)
plt.xlabel('Spanwise Position z [m]', fontsize=12)
plt.ylabel('Reserve Factor (Strength / Load)', fontsize=12)
plt.grid(True, which='both', linestyle='--')
plt.legend()

# Limit y-axis to focus on the critical area (near 1.0)
# Adjust the top limit (e.g., 5 or 10) depending on how high your margins go
plt.ylim(0, 10) 

plt.tight_layout()
plt.show()