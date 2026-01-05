import numpy as np
import matplotlib.pyplot as plt

#---Constants---------------------------------------------------------------
# Note: Python uses ** for exponents, not ^. 
# Alternatively, use scientific notation like 71e9
E = 71 * 10**9     # Pa (Young's Modulus)
v = 0.33           # Poisson's ratio
t = 0.005          # Skin thickness in meters 
L = 9.79          # Total Span of the wing in meters (EXAMPLE VALUE - CHANGE THIS)

#---Design Parameters Setup-------------------------------------------------
# 1. Create the Z coordinate array (500 data points from root to tip)
z = np.linspace(0, L, 500)

# 2. Create arrays for b and k_c that match the size of z
b = np.zeros_like(z)   # Creates an empty array of zeros
k_c = np.zeros_like(z) # Creates an empty array of zeros

#---USER INPUT SECTION------------------------------------------------------
# THIS IS WHERE YOU PUT YOUR VALUES
# --------------------------------------------------------------------------

# Define the spacing values for the two sections

w_wingbox = 2.8735 * 0.3 #width of the wingbox in meters 
#number of stringers for each design (half wing)
n_s_root1 = 4
n_s_tip1 = 2

n_s_root2 = 7
n_s_tip2  = 4

n_s_root3 = 9
n_s_tip3  = 5

b_value_root = (n_s_root1+1)/w_wingbox  # Stringer spacing for the first half (in meters)
b_value_tip  = (n_s_tip1+1)/w_wingbox # Stringer spacing for the second half (in meters)

# Define k_c (Buckling coefficient)
# If k_c is constant (4) everywhere:
k_c[:] = 4 
# OR if k_c also changes halfway, uncomment the lines below:
# k_c[z < L/2] = 4.0
# k_c[z >= L/2] = 3.6

# --------------------------------------------------------------------------
# Applying the logic: b changes halfway through the wing
# We use a mask: "z < L/2" finds all indices where we are in the first half
midpoint = L / 2

b[z < midpoint]  = b_value_root  # First half gets root spacing
b[z >= midpoint] = b_value_tip   # Second half gets tip spacing

#---Calculation-------------------------------------------------------------
# We use np.pi to allow array calculation
# The formula is applied element-wise to the entire array at once
Sigma_cr_skin = (k_c * (np.pi ** 2) * E) / (12 * (1 - v ** 2)) * (t / b) ** 2

#---Plotting----------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(z, Sigma_cr_skin / 10**6, label='Skin Critical Buckling Stress', color='blue') # Converted to MPa for readability

plt.title('Skin Buckling Stress distribution along Span')
plt.xlabel('Spanwise position z [m]')
plt.ylabel('Critical Buckling Stress [MPa]')
plt.grid(True, which='both', linestyle='--')
plt.legend()

# Add a vertical line to show where the switch happens
plt.axvline(x=midpoint, color='red', linestyle=':', label='Section Change')

plt.show()