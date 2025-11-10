import numpy as np
import math as m

import sys
import os

# Get the absolute path to the parent directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Add the parent directory to Python's module search path
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Now import your module
import variables.fixed_values

C_d0 = variables.fixed_values.C_d0
AR = variables.fixed_values.AR
e = variables.fixed_values.oswald_efficiency
B = variables.fixed_values.B
p_ISO = variables.fixed_values.p_ISO
T_ISO = variables.fixed_values.T_ISO
rho_ISO = variables.fixed_values.rho_ISO

## needs to be specified:
mass_fraction_climb = variables.fixed_values.mass_fraction_climb
climb_rate_requirement = variables.fixed_values.climb_rate_requirement

velocity_climb_rate = []
Mach_number = []
total_pressure = []
lapse_rate = []

wing_loading = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]


def climb_rate(wing_loading_input):
    global velocity_climb_rate, Mach_number, total_pressure, lapse_rate
    output = np.array([])
    
    velocity_climb_rate = velocity_climb_rate_func(wing_loading_input)
    Mach_number = Mach_number_func(velocity_climb_rate)
    total_pressure = total_pressure_func(Mach_number)
    lapse_rate = lapse_rate_func(total_pressure, Mach_number)

    for i in range(len(wing_loading_input)):
        output = np.append(output, (mass_fraction_climb/lapse_rate[i]) * (climb_rate_requirement/velocity_climb_rate[i] + (2 * m.sqrt(C_d0/(m.pi * AR * e)))))
    return output

def velocity_climb_rate_func(wing_loading_input):
    output = np.array([])
    for value in wing_loading_input:
        output = np.append(output, m.sqrt(value * 2/(rho_ISO*0.728)))
    return output
        
def Mach_number_func(velocity_climb_rate_output):
    output = np.array([])
    for value in velocity_climb_rate_output:
        output = np.append(output, value/m.sqrt(T_ISO*1.4*287))
    return output

def total_pressure_func(Mach_number_output):
    output = np.array([])
    for i, value in enumerate(Mach_number_output):
        output = np.append(output, 1 + 0.2 * Mach_number_output[i]**2)
    return output

def lapse_rate_func(total_pressure_output, Mach_number_output):
    output = np.array([])
    for i in range(len(total_pressure_output)):
        output = np.append(output, total_pressure_output[i] * (1 - (0.43 + 0.014 * B) * m.sqrt(Mach_number_output[i])))
    return output