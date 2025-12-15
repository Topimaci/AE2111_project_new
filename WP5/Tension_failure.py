import numpy as np
import matplotlib.pyplot as plt

### Values of moment along span
M_vals = np.load("M_vals.npy")
x_grid = np.load("X_grid.npy")
I_xx = np.load("I_xx.npy")
h_fs = np.load("h_front_spar.npy")
h_rs = np.load("h_rear_spar.npy")

#distance from the top spar to the neutral axis
x_c = (h_rs ** 2 + h_fs ** 2 + h_fs * h_rs) / (3 * (h_rs + h_fs))

#for POSITIVE LOAD CASES
#distance from the neutral axis to the lower left point of cross section
y_norm_stress_front = h_fs - x_c

stress = M_vals*y_norm_stress_front/I_xx


#print(y_norm_stress_front)
#ultimate stress = 510MPA
#yield stress = 450 MPA
stress_critical = 450000000   ###has to be discussed what we define as critical

stress_critical_array= np.full_like(x_grid, stress_critical)


plt.figure(figsize=(8,5))
plt.plot(x_grid, stress, label='Stress', color='orange')
plt.xlabel('Spanwise Location y [m]')
plt.ylabel('Stress')
plt.title('Stress along Span')
plt.grid(True)
plt.legend()
plt.show()



cutoff_value = 400
cutoff_stress = stress[cutoff_value]
range_value = 500 - cutoff_value
for i in range(range_value):
    print(i)
    stress[cutoff_value+i] = cutoff_stress

print(stress[cutoff_value])
print(stress[499])

margin_of_safety = stress_critical_array/stress

plt.figure(figsize=(8,5))
plt.plot(x_grid, margin_of_safety, label='Margin of saftey', color='orange')
plt.xlabel('Spanwise Location y [m]')
plt.ylabel('Margin of saftey')
plt.title('Margin of safety')
plt.grid(True)
plt.legend()
plt.show()



