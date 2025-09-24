import Planform_DESIGN1 as Pl
import fixed_values as fv
import math

start_pos_chord = 4
end_pos_chord = 6

sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE = Pl.calculate_geometric_parameters_wing(30.46, fv.AR, 0.68, 24.5)

sweep_LE = math.radians(sweep_LE)

# trailing edge sweep
sweep_TE = math.atan(math.tan(sweep_LE) + (c_tip - c_root) / (b / 2) )

sweep_LE = math.radians(sweep_LE)

x_1 = end_pos_chord - start_pos_chord
area_1 = (x_1) ** 2 * math.tan(sweep_LE) / 2 
area_2 = (x_1) * (c_root - x_1 * math.tan(sweep_LE))
area_3 = x_1 ** 2 * math.tan(sweep_TE) / 2 

total_area = area_1 + area_2 + area_3

print(total_area)

