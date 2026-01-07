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

print("front spar", h_fs[0])
print("rear spar", h_rs[0])

#for POSITIVE LOAD CASES
#distance from the neutral axis to the lower left point of cross section
y_norm_stress_front = h_fs - x_c

print("y distance", y_norm_stress_front[0])
print("I_xx", I_xx[0])
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

print(y_norm_stress_front[0])
print("Root moment", M_vals[0])
print("Root inertia",I_xx[0])

cutoff_value = 350
cutoff_stress = stress[cutoff_value]
range_value = 500 - cutoff_value
for i in range(range_value):
    stress[cutoff_value+i] = cutoff_stress

 

margin_of_safety = stress_critical_array/stress

#np.save("Design1_Tension_MoS", margin_of_safety)
#np.save("Design2_Tension_MoS", margin_of_safety)
#np.save("Design3_Tension_MoS", margin_of_safety)
#np.save("Design4_Tension_MoS", margin_of_safety)
#np.save("Design5_Tension_MoS", margin_of_safety)

#MoS_1 = np.load("Design1_Tension_MoS.npy")
#MoS_2 = np.load("Design2_Tension_MoS.npy")
#MoS_3 = np.load("Design3_Tension_MoS.npy")
#MoS_4 = np.load("Design4_Tension_MoS.npy")
#MoS_5 = np.load("Design5_Tension_MoS.npy")


print("Margin of safety root", margin_of_safety[0])
print("Margin of safety at cutoff value", margin_of_safety[400])
print("distance at which MoS is getting cut off", x_grid[400])
plt.figure(figsize=(8,5))


#plt.plot(x_grid, MoS_1, label='MoS Design 1', color='red')
#plt.plot(x_grid, MoS_2, label='MoS Design 2', color='green')
plt.plot(x_grid, margin_of_safety, label='MoS Design 3', color='blue')
#plt.plot(x_grid, MoS_3, label='Margin of safety', color='blue')
#plt.plot(x_grid, MoS_4, label='Margin of safety', color='blue')
#plt.plot(x_grid, MoS_5, label='Margin of safety', color='blue')




# horizontal dotted line at y = 1
plt.axhline(y=1.25, color='red', linestyle='--', label='Safety threshold')

plt.xlabel('Spanwise Location y [m]')
plt.ylabel('Margin of safety')
plt.title('Margin of safety')
plt.grid(True)
plt.legend()

#starts the graph from, (0,0) otherwise made graph look weird
plt.ylim(bottom=0)  

plt.show()


