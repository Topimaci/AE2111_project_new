import Planform_DESIGN1 as Pl
import fixed_values as fv
import math

sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE = Pl.calculate_geometric_parameters_wing(30.46, fv.AR, 0.68)

sweep_LE = math.radians(sweep_LE)

sweep_TE = math.atan(math.tan(sweep_LE) - (2 * (1 - taper) * c_root / b ))

sweep_TE = math.degrees(sweep_TE)

print(sweep_TE)
# trailing edge sweep

