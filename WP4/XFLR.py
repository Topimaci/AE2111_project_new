import numpy as np

# Load all numeric data after 20 header lines
data = np.genfromtxt("angleofattack0.txt", skip_header=21, skip_footer=60, invalid_raise=False)

# Assign only the columns you care about to separate variables
y_span = data[:, 0]
chord  = data[:, 1]
Ai     = data[:, 2]
Cl     = data[:, 3]
ICd    = data[:, 5]
Cm     = data[:, 7]

# Quick check
print("y positions:", y_span[:])
print("Chord lengths:", chord[:])
print("Lift coefficients:", Cl[:])
print("Induced drag:", ICd[:])
print("Quarter-chord moment:", Cm[:])

