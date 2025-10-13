## model, thrust[N], mass[kg], length, diameter, CG
engines = [
    ["Honeywell HTF7500E", 33486, 687, 2.282, 1.116],
    ["Honeywell HTF7350 turbofans", 32574, 696, 2.281, 1.156],
    ["Honeywell HTF7700L", 34100, 696, 2.282, 1.156],
    ["Pratt & Whitney Canada PW305A", 23200, 452.5, 1.892, 1.147],
    ["Honeywell HTF7250G", 33900, 696, 2.282, 1,156]
]

engines_sorted = []

## Sort the engines based on thrust
engines_sorted = sorted(engines, key=lambda x: x[1])

## Find the smallest one to surpass the thrust required (full thrust)
def engine_required(thrust_required):
    for engine in engines_sorted:
        if engine[1] >= thrust_required/2:
            return engine ## returns a list of all specific engine values

print(engines_sorted)