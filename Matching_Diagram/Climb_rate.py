import numpy as np
import math as m

C_d0 = 0.059
AR = 10
e = 0.8
B = 4
p_ISO = 101325
T_ISO = 273.15
rho_ISO = 1.225225

## needs to be specified:
mass_fraction_climb = 0.97
climb_rate_requirement = 15  # m/s

velocity_climb_rate = []
Mach_number = []
total_pressure = []
lapse_rate = []

def Climb_rate(wing_loading_input):
    output = []
    velocity_climb_rate = velocity_climb_rate(wing_loading_input)
    Mach_number = Mach_number(velocity_climb_rate)
    total_pressure = total_pressure(Mach_number)
    lapse_rate = lapse_rate(total_pressure, Mach_number)

    for i in range(len(wing_loading_input)):
        output = np.append((mass_fraction_climb/lapse_rate[i]) * (climb_rate_requirement/velocity_climb_rate[i] + (2 * m.sqrt(C_d0/(m.pi() * AR * e)))))

def velocity_climb_rate(wing_loading_input):
    output = []
    for value in wing_loading_input:
        output = np.append(output, m.sqrt(value**2/(rho_ISO*0.728)))
    return output
        
def Mach_number(velocity_climb_rate_output):
    output = []
    for value in velocity_climb_rate_output:
        output = np.append(output, value/m.sqrt(T_ISO*1.4*287))
    return output

def total_pressure(Mach_number_output):
    output = []
    for i, value in enumerate(Mach_number_output):
        output = np.append(output, 1 + 0.2 * Mach_number_output[i]**2)
    return output

def lapse_rate(total_pressure_output, Mach_number_output):
    output = []
    for i in range(len(total_pressure_output)):
        output = np.append(output, total_pressure_output[i] * (1 - (0.43 + 0.014 * B) * m.sqrt(Mach_number_output[i])))
    return output