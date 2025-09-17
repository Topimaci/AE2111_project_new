rho_at_airport = 1.225225
bypass_ratio = 4
V_approach = 66
C_l_max_landing = 2.3

def Minimum_speed_calculation_function(rho_at_airport, bypass_ratio, V_approach, C_l_max_landing):
    Wing_loading = 1/bypass_ratio * (rho_at_airport/2) * (V_approach/1.23)**2 * C_l_max_landing
    return Wing_loading

def Minimum_speed_function(Loading_array_input):
    Loading_array_output = []
    for value in Loading_array_input:
        Loading_array_output.append(Minimum_speed_calculation_function(rho_at_airport, bypass_ratio, V_approach, C_l_max_landing))
    
    return Loading_array_output

print(Minimum_speed_function([1, 2, 3, 4, 5]))
