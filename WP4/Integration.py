from TorqueDist import compute_case
from XFLR import y_span0, chord0, Cl0, ICd0, Cm0, ICd10, Cm10, Cl10


V_inf = 58
rho = 1.225
aoa_deg = 0.0


res0 = compute_case(y_span0, chord0, Cl0, Cl10, aoa_deg, ICd0 , ICd10, Cm0, Cm10, V_inf, rho)

x_grid  = res0["x_grid"]
T_total = res0["T_total"]

print(x_grid[:5], T_total[:5])

