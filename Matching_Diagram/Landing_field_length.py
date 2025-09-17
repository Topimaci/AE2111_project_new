mass_fraction1 = 0.7
landing_field_length1 = 950
density_land1 = 1.165
C_L_max_land1 = 2.3


def landing_field_length(mass_fraction, landing_field_length, density_land, C_L_max_land):
    Wing_loading = 1 / mass_fraction * landing_field_length * 1/0.9 * density_land * C_L_max_land
    return Wing_loading
