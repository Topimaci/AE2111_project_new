def Minimum_speed(rho_at_airport, bypass_ratio, V_approach, C_l_max_landing):
    Wing_loading = 1/bypass_ratio * (rho_at_airport/2) * (V_approach/1.23)**2 * C_l_max_landing
    return Wing_loading
