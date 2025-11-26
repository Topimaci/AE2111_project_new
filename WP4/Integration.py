from TorqueDist import compute_case
from XFLR import y_span0, chord0, Cl0, ICd0, Cm0

V_inf = 53
rho = 1.225

res0 = compute_case(y_span0, chord0, Cl0, ICd0, Cm0, aoa_deg_case=0.0, V_inf=V_inf, rho=rho)

x_grid  = res0["x_grid"]
T_total = res0["T_total"]

print(x_grid[:5], T_total[:5])
