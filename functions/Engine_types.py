

## Find the smallest one to surpass the thrust required (full thrust)
def engine_required(thrust_required):

    # Format: ["Engine name", T_per_engine_N, Weight_per_engine_kg, Wtot_plus_nacelle, Nacelle_diameter_m, 
    #           Nacelle_length_m, CG_from_front_m, B_ratio, Specific_fuel_consumption_cruise_kg/N/s]
    engines = [
    ["Honeywell HTF7500E",            33486, 687,   916,     1.116, 2.282, 0.9128, 4.2, 18.2],
    ["Honeywell HTF7350 turbofans",   32574, 696,   928,     1.156, 2.281, 0.9124, 4.2, 18.2],
    ["Honeywell HTF7700L",            34100, 696,   928,     1.156, 2.282, 0.9128, 4.2, 18.2],
    ["Pratt & Whitney Canada PW305A", 23200, 452.5, 603.333, 1.147, 2.06,  0.824,  4.5, 19.4],
    ["Honeywell HTF7250G",            33900, 696,   928,     1.156, 2.281, 0.9124, 4.2, 18.2]
]   


    engines_sorted = []

## Sort the engines based on thrust
    engines_sorted = sorted(engines, key=lambda x: x[1])


    for engine in engines_sorted:
            if engine[1] >= thrust_required/2:
                return engine  # found the smallest that meets requirement

        # If none meet the thrust requirement
    return ["00"]