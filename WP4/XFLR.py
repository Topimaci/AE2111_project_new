import numpy as np
import warnings

# Load all numeric data after 20 header lines (suppress warnings about missing columns)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    data = np.genfromtxt("angleofattack0.txt", skip_header=21, skip_footer=60, invalid_raise=False)

# In the XFLR file there are many columns that aren't of interest, for this reason only those needed are extracted from the data matrix found above. Keeping the data matrix allows 
# extracting other values if they are ever needed. first value ":" means all values in column (all rows), while second means "column number" (column position in the XFLR txt file).
y_span = data[:, 0]
chord  = data[:, 1]
Ai     = data[:, 2]
Cl     = data[:, 3]
ICd    = data[:, 5]
Cm     = data[:, 7]

#check
print("y positions:", y_span[:])
print("Chord lengths:", chord[:])
print("Lift coefficients:", Cl[:])
print("Induced drag:", ICd[:])
print("Quarter-chord moment:", Cm[:])
