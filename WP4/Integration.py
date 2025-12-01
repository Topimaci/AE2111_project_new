from TorqueDist import compute_case
from XFLR import y_span0, chord0, Cl0, ICd0, Cm0, ICd10, Cm10, Cl10
import numpy as np

V_inf = 200.736
rho = 0.3662

res0 = compute_case(y_span0, chord0, Cl0, Cl10, aoa_deg, ICd0 , ICd10, Cm0, Cm10, V_inf, rho)

x_grid  = np.array(res0["x_grid"])
T_total = np.array(res0["T_total"])


types = [type(x) for x in x_grid]
