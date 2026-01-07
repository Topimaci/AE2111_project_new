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
stress_applied = np.abs(M_vals * y_norm_stress_front / I_xx)

# 4. Apply Cutoff (Fixing the tip singularity)
cutoff_idx = 400
if len(stress_applied) > cutoff_idx:
    cutoff_stress = stress_applied[cutoff_idx]
    stress_applied[cutoff_idx:] = cutoff_stress

# ==============================================================================
# PART 2 & 3: CALCULATE CRITICAL STRESS & PLOT FOR ALL 5 DESIGNS
# ==============================================================================

# Constants
E = 71 * 10**9     # Pa
v = 0.33           # Poisson's ratio
z = x_grid
L = z[-1]          # Total span
w_wingbox = 2.8735 * 0.3  # width of wingbox [m]

# --- DEFINE ALL DESIGNS ---
designs = [
    # DESIGN 1 (Simple Split)
    {
        "label": "Design 1", 
        "type": "simple",
        "t": 0.002, 
        "n_root": 4, 
        "n_tip": 2, 
        "color": "blue"
    },
    # DESIGN 2 (Simple Split)
    {
        "label": "Design 2", 
        "type": "simple",
        "t": 0.003, 
        "n_root": 7, 
        "n_tip": 4, 
        "color": "orange"
    },
    # DESIGN 3 (Simple Split)
    {
        "label": "Design 3", 
        "type": "simple",
        "t": 0.003, 
        "n_root": 9, 
        "n_tip": 5, 
        "color": "green"
    },
    # DESIGN 4 (Updated: Spar/Skin Heavy)
    {
        "label": "Design 4 (Thick Skin)",
        "type": "complex",
        "t": 0.005,
        "y_breaks": np.array([0, 3, 4.89, 7]),
        "n_top": np.array([8, 7, 5, 4]), # Updated counts
        "color": "red"
    },
    # DESIGN 5 (Updated: Stringer Heavy)
    {
        "label": "Design 5 (Many Stringers)",
        "type": "complex",
        "t": 0.004,
        "y_breaks": np.array([0, 3, 4.89, 7]),
        "n_top": np.array([9, 8, 6, 4]), # Updated counts
        "color": "purple"
    }
]

# Initialize Plot
plt.figure(figsize=(14, 8))

# Loop through each design to calculate and plot
for d in designs:
    t_current = d["t"]
    
    # --- LOGIC TO GENERATE 'b' ARRAY ---
    b = np.zeros_like(z)
    
    if d["type"] == "simple":
        # Logic: Split halfway
        n_root = d["n_root"]
        n_tip = d["n_tip"]
        
        b_value_root = w_wingbox / (n_root + 1)
        b_value_tip  = w_wingbox / (n_tip + 1)
        
        midpoint = L / 2
        b[z < midpoint]  = b_value_root
        b[z >= midpoint] = b_value_tip
        
    elif d["type"] == "complex":
        # Logic: Use y_breaks to determine number of stringers
        breaks = d["y_breaks"]
        counts = d["n_top"]
        
        # Find which interval each z value falls into
        indices = np.searchsorted(breaks, z, side='right') - 1
        
        # Handle boundaries
        indices[indices < 0] = 0
        indices[indices >= len(counts)] = len(counts) - 1
        
        # Get count for every point z
        current_n_counts = counts[indices]
        
        # Calculate b
        b = w_wingbox / (current_n_counts + 1)

    # --- CALCULATE MARGIN ---
    k_c = 4.0  # Assume simply supported

    # Critical Buckling Stress
    Sigma_cr_skin = (k_c * (np.pi ** 2) * E) / (12 * (1 - v ** 2)) * (t_current / b) ** 2

    # Reserve Factor
    with np.errstate(divide='ignore', invalid='ignore'):
        margin_of_safety = Sigma_cr_skin / stress_applied

    # Plot
    plt.plot(z, margin_of_safety, label=f'{d["label"]} (t={t_current*1000}mm)', color=d["color"], linewidth=2)


# --- FORMATTING THE GRAPH ---

# Safety Thresholds
plt.axhline(y=1.0, color='red', linestyle='--', linewidth=3, label='Failure Threshold (RF=1)')
plt.axhline(y=1.25, color='black', linestyle=':', linewidth=2, label='Safety Target (RF=1.25)')

plt.title('Skin Buckling Safety Margin Comparison (All 5 Designs)', fontsize=14)
plt.xlabel('Spanwise Position z [m]', fontsize=12)
plt.ylabel('Reserve Factor (Strength / Load)', fontsize=12)
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(loc='upper right')

# Axis limits
plt.ylim(0, 8) 

plt.tight_layout()
plt.show()