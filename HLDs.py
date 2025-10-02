import Planform_DESIGN1 as Pl
import fixed_values as fv
import math

start_pos_span_TE = 0.15
end_pos_span_TE = 0.8

start_pos_span_LE = 0.2
end_pos_span_LE = 0.8
c_f_c_TE = 0.3
c_f_c_LE = 0.1

sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE = Pl.calculate_geometric_parameters_wing(30.46, fv.AR, 0.68, 24.5)

sweep_LE = math.radians(sweep_LE)
half_b = b / 2 

# trailing edge reference area
sweep_TE = math.atan(math.tan(sweep_LE) + (c_tip - c_root) / (b / 2) )

length_TE = (end_pos_span_TE - start_pos_span_TE) * b / 2
area_1 = (length_TE) ** 2 * math.tan(sweep_LE) / 2 
area_2 = (length_TE) * (c_root - length_TE * math.tan(sweep_LE))
area_3 = length_TE ** 2 * math.tan(sweep_TE) / 2 

reference_area_TE = (area_1 + area_2 + area_3) * 2

# leading edge reference area
Datum = c_tip + half_b * math.tan(sweep_LE)
y_1 = start_pos_span_LE * half_b * math.tan(sweep_LE)
x_1 = end_pos_span_LE * half_b * math.tan(sweep_LE)
x_2 = end_pos_span_LE * half_b * math.tan(sweep_TE)
y_2 = start_pos_span_LE * half_b * math.tan(sweep_TE)

c_1 = Datum - x_1 - x_2 
c_2 = Datum - y_1 - y_2 

reference_area_LE = (c_1 + c_2) / 2 * (end_pos_span_LE - start_pos_span_LE) * half_b


# C_L_max calculation trailing edge fowler flap

c_avg_TE = c_root - (end_pos_span_TE - start_pos_span_TE) / 2 * b * math.tan(sweep_LE)
delc_cf_TE = 0.6
c_prime_TE = c_avg_TE + delc_cf_TE * c_f_c_TE * c_avg_TE
delta_c_lmax_TE = 1.3 * c_prime_TE / c_avg_TE
wing_area = (c_root + c_tip) / 2 * b 

DELTA_c_LMAX_TE = 0.9 * delta_c_lmax_TE * reference_area_TE / wing_area * math.cos(sweep_TE)


# C_L_max calculation leading edge slats

b_avg = (start_pos_span_LE + x_2 / 2) * half_b  # average position of leading edge HLD
c_avg_LE = c_root - math.tan(sweep_LE) * b_avg + b_avg * math.tan(sweep_TE)
delc_cf_LE = 0.2
c_prime_LE = c_avg_LE + delc_cf_LE * c_f_c_LE * c_avg_LE
delta_c_lmax_LE = 0.4 * c_prime_LE / c_avg_LE
DELTA_c_LMAX_LE = 0.9 * delta_c_lmax_LE * reference_area_LE / wing_area * math.cos(sweep_LE)

print(DELTA_c_LMAX_LE + DELTA_c_LMAX_TE)

print(reference_area_LE)