import numpy as np
import matplotlib.pyplot as plt

# ==============================================================================
# PART 1: LOAD DATA
# ==============================================================================
M_vals = np.load("M_vals.npy")
x_grid = np.load("X_grid.npy")
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
    y_dist = h_avg / 2.0
    
    # 3. Calculate Inertia Components
    # A. Spars (Vertical webs) - I = 1/12 * t * h^3 (for two spars)
    I_spars = (1/12 * t_spar * h_fs**3) + (1/12 * t_spar * h_rs**3)
    
    # B. Skins (Horizontal plates) - I = A * y^2 (Parallel Axis Theorem)
    I_skins = 2 * (w * t_skin) * y_dist**2
    
    # C. Stringers (Lumped areas) - I = n * A * y^2
    I_stringers_top = n_top_arr * A_str * y_dist**2
    I_stringers_bot = n_bot_arr * A_str * y_dist**2
    
    # Total I_xx
    I_xx_total = I_spars + I_skins + I_stringers_top + I_stringers_bot
    
    return I_xx_total, y_dist

# ==============================================================================
# PART 3: DEFINE DESIGNS
# ==============================================================================
designs = [
    # DESIGN 1
#    {
#        "label": "Design 1 (Light)", "type": "simple", "color": "blue",
 #       "t": 0.002, 
  #      "t_spar": 0.005, "A_str": 0.0001,
   #     "n_root": 4, "n_tip": 2,
    #    "n_bot_root": 4, "n_bot_tip": 2
    #},
    # DESIGN 2
    #{
     #   "label": "Design 2 (Med)", "type": "simple", "color": "orange",
      #  "t": 0.003, 
       # "t_spar": 0.008, "A_str": 0.00025,
        #"n_root": 7, "n_tip": 4,
        #"n_bot_root": 7, "n_bot_tip": 4
    #},
    # DESIGN 3
    #{
     #   "label": "Design 3 (Heavy)", "type": "simple", "color": "green",
      #  "t": 0.003, 
       # "t_spar": 0.005, "A_str": 0.0002,
        #"n_root": 9, "n_tip": 5,
        #"n_bot_root": 9, "n_bot_tip": 5
    #},
    # DESIGN 4.0
    {
        "label": "Design 4 (Thick Skin)", "type": "complex", "color": "red",
        "t": 0.005, 
        "t_spar": 0.010, "A_str": 0.0004,
        "y_breaks": np.array([0, 3, 4.89, 7]),
        "n_top": np.array([8, 7, 5, 4]), 
        "n_bot": np.array([5, 5, 4, 4]) 
    },
    # DESIGN 5
    {
        "label": "Design 5 (Many Stringers)", "type": "complex", "color": "purple",
        "t": 0.004, 
        "t_spar": 0.008, "A_str": 0.0003,
        "y_breaks": np.array([0, 3, 4.89, 7]),
        "n_top": np.array([9, 8, 6, 4]), 
        "n_bot": np.array([7, 7, 4, 3]) 
    }
]

# ==============================================================================
# PART 4: MAIN CALCULATION & PLOTTING LOOP
# ==============================================================================
plt.figure(figsize=(14, 8))
cutoff_idx = 400 

for d in designs:
    # --- A. DETERMINE STRINGER DISTRIBUTION (n_top and n_bot) ---
    n_top_arr = np.zeros_like(z)
    n_bot_arr = np.zeros_like(z)
    
    if d["type"] == "simple":
        midpoint = L / 2
        mask_root = z < midpoint
        mask_tip  = z >= midpoint
        
        n_top_arr[mask_root] = d["n_root"]
        n_top_arr[mask_tip]  = d["n_tip"]
        n_bot_arr[mask_root] = d["n_bot_root"]
        n_bot_arr[mask_tip]  = d["n_bot_tip"]
        
    elif d["type"] == "complex":
        indices = np.searchsorted(d["y_breaks"], z, side='right') - 1
        indices[indices < 0] = 0
        indices[indices >= len(d["n_top"])] = len(d["n_top"]) - 1
        
        n_top_arr = d["n_top"][indices]
        n_bot_arr = d["n_bot"][indices]

    # --- B. CALCULATE I_XX & APPLIED STRESS ---
    # Note: I_xx depends on BOTH top and bottom stringers
    I_xx_custom, y_max = calculate_structure(z, h_fs, h_rs, w_wingbox, d["t"], d["t_spar"], n_top_arr, n_bot_arr, d["A_str"])
    
    stress_applied = np.abs(M_vals * y_max / I_xx_custom)
    
    # Fix tip singularity
    if len(stress_applied) > cutoff_idx:
        stress_applied[cutoff_idx:] = stress_applied[cutoff_idx]

    # --- C. TOP SKIN ANALYSIS (Standard for all designs) ---
    b_top = w_wingbox / (n_top_arr + 1)
    k_c = 4.0
    
    Sigma_cr_top = (k_c * (np.pi ** 2) * E) / (12 * (1 - v ** 2)) * (d["t"] / b_top) ** 2
    
    with np.errstate(divide='ignore', invalid='ignore'):
        margin_top = Sigma_cr_top / stress_applied

    # Plot TOP Margin
    plt.plot(z, margin_top, label=f'{d["label"]} (Top)', color=d["color"], linewidth=2)

    # --- D. BOTTOM SKIN ANALYSIS (Only for Design 4 and 5) ---
    if d["label"].startswith("Design 4") or d["label"].startswith("Design 5"):
        
        # Calculate spacing for bottom skin
        b_bot = w_wingbox / (n_bot_arr + 1)
        
        # Calculate Critical Stress for Bottom
        Sigma_cr_bot = (k_c * (np.pi ** 2) * E) / (12 * (1 - v ** 2)) * (d["t"] / b_bot) ** 2
        
        # Calculate Margin for Bottom
        # Assumption: Checking bottom skin buckling against the same compressive load magnitude
        with np.errstate(divide='ignore', invalid='ignore'):
            margin_bot = Sigma_cr_bot / stress_applied
            
        # Plot BOTTOM Margin (Dashed Line)
        plt.plot(z, margin_bot, label=f'{d["label"]} (Bottom)', color=d["color"], linestyle='--', linewidth=2)


# --- FORMATTING ---
plt.axhline(y=1.0, color='red', linestyle='--', linewidth=3, label='Failure Threshold (RF=1.0)')
plt.axhline(y=1.25, color='black', linestyle=':', linewidth=2, label='Safety Target (RF=1.25)')

plt.title('Skin Buckling Safety Margin (Top vs Bottom for D4/D5)', fontsize=16)
plt.xlabel('Spanwise Position z [m]', fontsize=14)
plt.ylabel('Safety Factor', fontsize=14)
plt.grid(True, which='both', linestyle='--', alpha=0.7)
plt.legend(loc='upper right', fontsize=10, ncol=2) # Two columns for legend to save space

plt.ylim(0, 10) 
plt.tight_layout()
plt.show()