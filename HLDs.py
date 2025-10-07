import Planform_DESIGN1 as Pl
import fixed_values as fv
import dynamic_variables as dv
import math

start_pos_span_TE = 0.15
end_pos_span_TE = 0.8

start_pos_span_LE = 0.2
end_pos_span_LE = 0.8
c_f_c_TE = 0.35
c_f_c_LE = 0.1
radius_fuselage = 1

sweep, taper, b, c_root, c_tip, c_MAC, dihedral, sweep_LE = Pl.calculate_geometric_parameters_wing(dv.S_w, fv.AR, 0.68, 24.5)


# General parameters
sweep_LE = math.radians(sweep_LE)
half_b = b / 2 
Datum = c_tip + half_b * math.tan(sweep_LE)

# trailing edge reference area
sweep_TE = math.atan(math.tan(sweep_LE) + (c_tip - c_root) / (b / 2) )
z_1 = radius_fuselage * math.tan(sweep_LE)
z_2 = (half_b - radius_fuselage) * math.tan(sweep_TE)
q_1 = end_pos_span_TE * half_b * math.tan(sweep_LE)
q_2 = (half_b - end_pos_span_TE * half_b) * math.tan(sweep_TE)

t_1 = Datum - q_1 - q_2
t_2 = Datum - z_1 - z_2

reference_area_TE = ((t_1 + t_2) / 2 * (end_pos_span_TE * half_b - radius_fuselage)) * 2


# leading edge reference area
y_1 = end_pos_span_LE * half_b * math.tan(sweep_LE)
y_2 = (half_b - end_pos_span_LE * half_b) * math.tan(sweep_TE)
x_1 = start_pos_span_LE * half_b * math.tan(sweep_LE)
x_2 = end_pos_span_LE * half_b * math.tan(sweep_TE)

c_1 = Datum - x_1 - x_2 
c_2 = Datum - y_1 - y_2 

reference_area_LE = ((c_1 + c_2) / 2 * (end_pos_span_LE - start_pos_span_LE) * half_b) * 2


# C_L_max calculation trailing edge fowler flap

c_avg_TE = Datum - c_tip - (end_pos_span_TE - start_pos_span_TE) / 2 * b * math.tan(sweep_LE)
delc_cf_TE = 0.6
c_prime_TE = c_avg_TE + delc_cf_TE * c_f_c_TE * c_avg_TE
delta_c_lmax_TE = 1.3 * c_prime_TE / c_avg_TE
wing_area = (c_root + c_tip) / 2 * b 
print(wing_area)

DELTA_c_LMAX_TE = 0.9 * delta_c_lmax_TE * reference_area_TE / wing_area * math.cos(sweep_TE)


# C_L_max calculation leading edge slats

b_avg = (start_pos_span_LE + x_2 / 2) * half_b  # average position of leading edge HLD
c_avg_LE = Datum - c_tip - math.tan(sweep_LE) * b_avg + b_avg * math.tan(sweep_TE)
delc_cf_LE = 0.2
c_prime_LE = c_avg_LE + delc_cf_LE * c_f_c_LE * c_avg_LE
delta_c_lmax_LE = 0.4 * c_prime_LE / c_avg_LE
DELTA_c_LMAX_LE = 0.9 * delta_c_lmax_LE * reference_area_LE / wing_area * math.cos(sweep_LE)

print(DELTA_c_LMAX_LE + DELTA_c_LMAX_TE)
print(c_avg_LE, c_avg_TE, b_avg)