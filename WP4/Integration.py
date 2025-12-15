from TorqueDist import compute_case, critical_alpha
from XFLR import y_span0, chord0, Cl0, ICd0, Cm0, ICd10, Cm10, Cl10
import numpy as np
import conditions as c

V_inf = c.velocity
rho = c.density
aoa_deg = critical_alpha(rho, V_inf, 38.379, c.weight, c.load_factor, c.landing, c.takeoff)
aoa_deg_abs = critical_alpha(rho, V_inf, 38.379, c.weight, abs(c.load_factor), c.landing, c.takeoff)

res0 = compute_case(y_span0, chord0, Cl0, Cl10, aoa_deg, ICd0 , ICd10, Cm0, Cm10, V_inf, rho)

x_grid  = np.array(res0["x_grid"])
T_total = np.array(res0["T_total"])


types = [type(x) for x in x_grid]

np.save("X_grid", x_grid)