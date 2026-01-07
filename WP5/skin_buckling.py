import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# PART 1: LOAD DATA
# ==============================================================================
# Ensure these files are in the same directory
M_vals = np.load("M_vals.npy")
x_grid = np.load("X_grid.npy")  # Master "z" coordinate
h_fs = np.load("h_front_spar.npy")
h_rs = np.load("h_rear_spar.npy")

# Constants
E = 71 * 10**9     # Pa
v = 0.33           # Poisson's ratio
z = x_grid
L = z[-1]          # Total span
w_wingbox = 2.8735 * 0.3  # width of wingbox [m]

# ==============================================================================
# PART 2: HELPER FUNCTION TO CALCULATE I_XX
# ==============================================================================
def calculate_structure(z_array, h_fs, h_rs, w, t_skin, t_spar, n_top_arr, n_bot_arr, A_str):
    """
    Calculates I_xx and y_max (distance to furthest fiber) for the wingbox
    at every point z.
    Approximation: Idealized Box with lumped stringers.
    """
    # 1. Average Height of the box at each station
    h_avg = (h_fs + h_rs) / 2.0
    
    # 2. Vertical distance from Neutral Axis (assuming symmetric-ish box)
    # y = h/2
    y_dist = h_avg / 2.0
    
    # 3. Calculate Inertia Components
    # A. Spars (Vertical webs) - I = 1/12 * t * h^3 (for two spars)
    I_spars = (1/12 * t_spar * h_fs**3) + (1/12 * t_spar * h_rs**3)
    
    # B. Skins (Horizontal plates) - I = A * y^2 (Parallel Axis Theorem)
    # Area of one skin = w * t_skin
    I_skins = 2 * (w * t_skin) * y_dist**2
    
    # C. Stringers (Lumped areas) - I = n * A * y^2
    I_stringers_top = n_top_arr * A_str * y_dist**2
    I_stringers_bot = n_bot_arr * A_str * y_dist**2
    
    # Total I_xx
    I_xx_total = I_spars + I_skins + I_stringers_top + I_stringers_bot
    
    # Return I_xx and the 'y' distance for stress calc
    return I_xx_total, y_dist

# ==============================================================================
# PART 3: DEFINE DESIGNS (UPDATED WITH YOUR VALUES)
# ==============================================================================
designs = [
    # DESIGN 1 (From Screenshot)
   # {
   #     "label": "Design 1 (Light)", "type": "simple", "color": "blue",
   #     "t": 0.002, 
   #     "t_spar": 0.005, "A_str": 0.0001,
   #     "n_root": 4, "n_tip": 2,
   #     "n_bot_root": 4, "n_bot_tip": 2
   # },
    # DESIGN 2 (From Screenshot)
   # {
   #     "label": "Design 2 (Med)", "type": "simple", "color": "orange",
   #     "t": 0.003, 
   #     "t_spar": 0.008, "A_str": 0.00025,
   #     "n_root": 7, "n_tip": 4,
   #     "n_bot_root": 7, "n_bot_tip": 4
   # },
    # DESIGN 3 (From Screenshot)
   # {
   #     "label": "Design 3 (Heavy)", "type": "simple", "color": "green",
   #     "t": 0.003, 
   #     "t_spar": 0.005, "A_str": 0.0002,
   #     "n_root": 9, "n_tip": 5,
   #     "n_bot_root": 9, "n_bot_tip": 5
   # },
    # DESIGN 4 (Updated Text - Spar/Skin Heavy)
    {
        "label": "Design 4 (Thick Skin)", "type": "complex", "color": "red",
        "t": 0.005, 
        "t_spar": 0.010, "A_str": 0.0004,
        "y_breaks": np.array([0, 3, 4.89, 7]),
        "n_top": np.array([8, 7, 5, 4]), 
        "n_bot": np.array([5, 5, 4, 4]) 
    },
    # DESIGN 5 (Updated Text - Stringer Heavy)
    {
        "label": "Design 5 (Many Stringers)", "type": "complex", "color": "purple",
        "t": 0.004, 
        "t_spar": 0.008, "A_str": 0.00035,
        "y_breaks": np.array([0, 3, 4.89, 7]),
        "n_top": np.array([9, 8, 6, 4]), 
        "n_bot": np.array([6, 6, 3, 3]) 
    }
]

# ==============================================================================
# PART 4: MAIN CALCULATION & PLOTTING LOOP
# ==============================================================================
plt.figure(figsize=(14, 8))

# Cutoff index to avoid tip singularity issues in plotting
cutoff_idx = 400 

for d in designs:
    # --- A. DETERMINE STRINGER DISTRIBUTION (n_top and n_bot) ---
    n_top_arr = np.zeros_like(z)
    n_bot_arr = np.zeros_like(z)
    b_spacing = np.zeros_like(z)
    
    if d["type"] == "simple":
        midpoint = L / 2
        mask_root = z < midpoint
        mask_tip  = z >= midpoint
        
        # Top Stringers
        n_top_arr[mask_root] = d["n_root"]
        n_top_arr[mask_tip]  = d["n_tip"]
        # Bottom Stringers
        n_bot_arr[mask_root] = d["n_bot_root"]
        n_bot_arr[mask_tip]  = d["n_bot_tip"]
        
    elif d["type"] == "complex":
        # Find indices for intervals based on y_breaks
        indices = np.searchsorted(d["y_breaks"], z, side='right') - 1
        indices[indices < 0] = 0
        indices[indices >= len(d["n_top"])] = len(d["n_top"]) - 1
        
        n_top_arr = d["n_top"][indices]
        n_bot_arr = d["n_bot"][indices]

    # Calculate Spacing 'b' based on Top Stringers (Critical for buckling)
    # Note: Adding 1 to stringer count gives number of bays
    b_spacing = w_wingbox / (n_top_arr + 1)
    
    # --- B. CALCULATE I_XX & APPLIED STRESS ---
    # Calculates specific Inertia for this design
    I_xx_custom, y_max = calculate_structure(z, h_fs, h_rs, w_wingbox, d["t"], d["t_spar"], n_top_arr, n_bot_arr, d["A_str"])
    
    # Calculate Stress: M * y / I
    stress_applied = np.abs(M_vals * y_max / I_xx_custom)
    
    # Fix tip singularity for cleaner plot
    if len(stress_applied) > cutoff_idx:
        stress_applied[cutoff_idx:] = stress_applied[cutoff_idx]

    # --- C. CALCULATE STRENGTH & MARGIN ---
    k_c = 4.0 # Simply supported assumption
    
    # Buckling Strength
    Sigma_cr_skin = (k_c * (np.pi ** 2) * E) / (12 * (1 - v ** 2)) * (d["t"] / b_spacing) ** 2
    
    # Reserve Factor (Margin of Safety)
    with np.errstate(divide='ignore', invalid='ignore'):
        margin_of_safety = Sigma_cr_skin / stress_applied

    # --- D. PLOT ---
    # Plotting up to the cutoff or end of span
    plt.plot(z, margin_of_safety, label=f'{d["label"]}', color=d["color"], linewidth=2)

# --- FORMATTING ---
plt.axhline(y=1.0, color='red', linestyle='--', linewidth=3, label='Failure Threshold (RF=1.0)')
plt.axhline(y=1.25, color='black', linestyle=':', linewidth=2, label='Safety Target (RF=1.25)')

plt.title('Skin Buckling Safety Margin Comparison (Dynamic I_xx)', fontsize=16)
plt.xlabel('Spanwise Position z [m]', fontsize=14)
plt.ylabel('Safety Factor', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(loc='upper right', fontsize=12)

# Adjust y-limits to focus on the critical region near RF=1
plt.ylim(0, 10) 

plt.tight_layout()
plt.show()