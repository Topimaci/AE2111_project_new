
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Shared data (same for all designs)
# -----------------------------
M_vals = np.load("M_vals.npy")
x_grid = np.load("X_grid.npy")

# -----------------------------
# Critical stress (compression)
# -----------------------------
stress_critical = 450000000  # Pa
stress_critical_array = np.full_like(x_grid, stress_critical)

# -----------------------------
# Cutoff parameters
# -----------------------------
cutoff_value = 350
range_value = 500 - cutoff_value

# -----------------------------
# first graph Designs 1, 2, 3
# -----------------------------
designs_123 = ["D1", "D2", "D3"]

plt.figure(figsize=(8,5))

for d in designs_123:

    I_xx = np.load(f"I_xx_{d}.npy")
    h_fs = np.load(f"h_front_spar_{d}.npy")
    h_rs = np.load(f"h_rear_spar_{d}.npy")

    x_c = (h_rs**2 + h_fs**2 + h_fs*h_rs) / (3*(h_rs + h_fs))
    stress = M_vals * x_c / I_xx
    
    # Apply cutoff
    cutoff_stress = stress[cutoff_value]
    for i in range(range_value):
        stress[cutoff_value + i] = cutoff_stress

    margin_of_safety = stress_critical_array / stress

    plt.plot(x_grid, margin_of_safety, label=f"Design {d}")

plt.axhline(y=1, color='red', linestyle='--', label='Safety threshold')
plt.xlabel('Spanwise Location y [m]')
plt.ylabel('Margin of Safety')
plt.title('Compression Margin of Safety - Designs 1, 2, 3')
plt.grid(True)
plt.legend()
plt.ylim(bottom=0)
plt.show()


# -----------------------------
#second graph, Designs 4, 5
# -----------------------------
designs_45 = ["D4", "D5"]

plt.figure(figsize=(8,5))

for d in designs_45:

    I_xx = np.load(f"I_xx_{d}.npy")
    h_fs = np.load(f"h_front_spar_{d}.npy")
    h_rs = np.load(f"h_rear_spar_{d}.npy")

    x_c = (h_rs**2 + h_fs**2 + h_fs*h_rs) / (3*(h_rs + h_fs))
    stress = M_vals * x_c / I_xx
    print(d) 
    print(I_xx)
    # Apply cutoff
    cutoff_stress = stress[cutoff_value]
    for i in range(range_value):
        stress[cutoff_value + i] = cutoff_stress

    margin_of_safety = stress_critical_array / stress

    plt.plot(x_grid, margin_of_safety, label=f"Design {d}")

plt.axhline(y=1, color='red', linestyle='--', label='Safety threshold')
plt.xlabel('Spanwise Location y [m]')
plt.ylabel('Margin of Safety')
plt.title('Compression Margin of Safety - Designs 4 & 5')
plt.grid(True)
plt.legend()
plt.ylim(bottom=0)
plt.show()
