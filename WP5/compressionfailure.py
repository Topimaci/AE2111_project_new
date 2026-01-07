import numpy as np
import matplotlib.pyplot as plt

# --- Load shared data ---
M_vals = np.load("M_vals.npy")
x_grid = np.load("X_grid.npy")

# --- Critical stress ---
stress_critical = 450_000_000  # Pa
stress_critical_array = np.full_like(x_grid, stress_critical)

# --- Cutoff parameters ---
cutoff_value = 350
range_value = 500 - cutoff_value

# --- Helper function to compute MoS ---
def compute_mos(I_file, h_fs_file, h_rs_file, save_prefix):
    I_xx = np.load(I_file)
    h_fs = np.load(h_fs_file)
    h_rs = np.load(h_rs_file)

    # Neutral axis to extreme fiber
    x_c = (h_rs ** 2 + h_fs ** 2 + h_fs * h_rs) / (3 * (h_rs + h_fs))

    # Bending stress
    stress = M_vals * (- x_c + h_fs) / I_xx

    # Apply cutoff
    max_index = min(cutoff_value + range_value, len(stress))
    stress[cutoff_value:max_index] = stress[cutoff_value]

    # Margin of safety
    mos = stress_critical_array / np.abs(stress)


    return mos

# --- Graph 1: D1, D2, D3 ---
plt.figure(figsize=(8,5))
for d in ["D1", "D2", "D3"]:
    mos = compute_mos(f"I_xx_{d}.npy", f"h_front_spar_{d}.npy", f"h_rear_spar_{d}.npy", save_prefix=d)
    plt.plot(x_grid[:len(mos)], mos, label=f"Design {d}")
plt.axhline(y=1.25, color='red', linestyle='--', label='Safety threshold')
plt.xlabel("Spanwise Location y [m]")
plt.ylabel("Margin of Safety")
plt.title("Compression Margin of Safety - Designs 1, 2, 3")
plt.grid(True)
plt.legend()
plt.ylim(bottom=0)
plt.show()

# --- Graph 2: D4, D5 ---
plt.figure(figsize=(8,5))
for d in ["D4", "D5"]:
    mos = compute_mos(f"I_xx_{d}.npy", f"h_front_spar_{d}.npy", f"h_rear_spar_{d}.npy", save_prefix=d)
    plt.plot(x_grid[:len(mos)], mos, label=f"Design {d}")
plt.axhline(y=1.25, color='red', linestyle='--', label='Safety threshold')
plt.xlabel("Spanwise Location y [m]")
plt.ylabel("Margin of Safety")
plt.title("Compression Margin of Safety - Designs 4 & 5")
plt.grid(True)
plt.legend()
plt.ylim(bottom=0)
plt.show()
