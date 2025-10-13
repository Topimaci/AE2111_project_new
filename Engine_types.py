## model, thrust, mass, length, diameter
engines = [
    ["eng 1", 1, 10, 100, 1000],
    ["eng 3", 3, 30, 300, 3000],
    ["eng 2", 2, 20, 200, 2000]
]

engines_sorted = []

## Sort the engines based on thrust
engines_sorted = sorted(engines, key=lambda x: x[1])

## Find the smallest one to surpass the thrust required
def engine_required(thrust_required):
    for engine in engines_sorted:
        if engine[1] >= thrust_required:
            return engine ## index
        

print(engine_required(0))