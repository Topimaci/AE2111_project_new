import numpy as np


#---Constants---------------------------------------------------------------

# --- Design 1 ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_1_S1 = 6.0
Stringers_1_S2 = 6.0
Stringers_1_S3 = 3.0
Stringers_1_S4 = 3.0

A_1 = 45.0 # Area stringer (crossection) in cm^2

# --- Design 2 ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_2_S1 = 4.0
Stringers_2_S2 = 4.0
Stringers_2_S3 = 2.0
Stringers_2_S4 = 2.0

A_2 = 65.0 # Area stringer (crossection) in cm^2

# --- Design 3 ---
# Number of stringers per section (0-3m, 3-4.87m, 4.87-7m, 5-9.8m)
Stringers_3_S1 = 4.0
Stringers_3_S2 = 3.0
Stringers_3_S3 = 2.0
Stringers_3_S4 = 2.0

A_3 = 65.0 # Area stringer (crossection) in cm^2

# --- Other constants ---
E = 71 * 10**9 # Young's modulus in Pa
k = 5.57 # Boundary condition constant (=4 due to the assumption that both ends are clamped)

