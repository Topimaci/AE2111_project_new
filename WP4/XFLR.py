import numpy as np
import warnings

# === AoA = 0° ===
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    data0 = np.genfromtxt(
        "WP4/data/angleofattack0.txt",
        skip_header=21,
        skip_footer=60,
        invalid_raise=False
    )

v = 200
M = v/343
y_span0 = data0[:, 0]
chord0  = data0[:, 1]
Ai0     = data0[:, 2]
Cl0     = data0[:, 3] *np.sqrt(1-M**2)
ICd0    = data0[:, 5]
Cm0     = data0[:, 7]

# (optional aliases if some other old code still expects these names)
y_span = y_span0
chord  = chord0
Ai     = Ai0
Cl     = Cl0
ICd    = ICd0
Cm     = Cm0

# === AoA = 10° ===
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=UserWarning)
    data10 = np.genfromtxt(
        "WP4/data/angleofattack10.txt",
        skip_header=21,
        skip_footer=60,
        invalid_raise=False
    )

y_span10 = data10[:, 0]
chord10  = data10[:, 1]
Ai10     = data10[:, 2]
Cl10     = data10[:, 3] * np.sqrt(1-M**2)
ICd10    = data10[:, 5]
Cm10     = data10[:, 7]

