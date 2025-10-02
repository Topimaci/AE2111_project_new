import Planform_DESIGN1 as Pl
import fixed_values as fv
import math

pos_span_TE = 0.4

start_pos_span_LE = 0.2
end_pos_span_LE = 0.55

sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE = Pl.calculate_geometric_parameters_wing(30.46, fv.AR, 0.68, 24.5)

sweep_LE = math.radians(sweep_LE)

# trailing edge reference area
sweep_TE = math.atan(math.tan(sweep_LE) + (c_tip - c_root) / (b / 2) )

x_1 = pos_span_TE * b / 2
area_1 = (x_1) ** 2 * math.tan(sweep_LE) / 2 
area_2 = (x_1) * (c_root - x_1 * math.tan(sweep_LE))
area_3 = x_1 ** 2 * math.tan(sweep_TE) / 2 

reference_area_TE = (area_1 + area_2 + area_3) * 2

# leading edge reference area
x_2 = end_pos_span_LE - start_pos_span_LE  
half_b = b / 2 
Area1 = x_2 * half_b * math.sin(sweep_LE) * x_2 * half_b / 2  # top triangle
Area2 = x_2 * half_b * (c_root - (x_2 + start_pos_span_LE) * half_b * math.sin(sweep_LE))  # big rectangle
Area3 = start_pos_span_LE * half_b * math.sin(sweep_TE) * x_2 * half_b  # small rectangle
Area4 = x_2 * half_b * math.sin(sweep_TE) * x_2 * half_b / 2  # bottom triangle

reference_area_LE = (Area1 + Area2 + Area3 + Area4) * 2

print(reference_area_LE, reference_area_TE)
