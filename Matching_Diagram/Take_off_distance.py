import numpy
import math

def find_V_2(wing_loading, density_take_off, C_L_take_off):
    V_2 = []
    for i in wing_loading: 
        x = math.sqrt(2 * i / density_take_off / C_L_take_off)
        V_2.append(x)
    return V_2

def find_Mach(V_2, temperature_take_off):
    Mach_take = []
    for j in V_2: 
        y = j / math.sqrt(1.4 * 287 * temperature_take_off)
        Mach_take.append(y)
    return Mach_take

def find_theta_delta(temperature_take_off, Mach_take):
    theta = []
    delta = []
    for k in Mach_take: 
        t = (temperature_take_off * (1 + 0.2 * k ** 2)) / 288.15
        d = (temperature_take_off * (1 + 0.2 * k ** 2) ** 3.5) / 101325
        theta.append(t)
        delta.append(d)
    return theta, delta

def find_alpha_t(delta, Mach_take, bypass_ratio):
    alpha_t = [l * (1 - (0.43 + 0.014 * bypass_ratio) * math.sqrt(m)) for l, m in zip(delta, Mach_take)]
    return alpha_t

def take_off_distance(alpha_t, wing_loading, take_off_field_length, density_take_off, oswald_efficiency, aspect_ratio):
    T_over_W = numpy.array()
    for wl, at in zip(wing_loading, alpha_t):
        q = 1 / at * (1.15 * math.sqrt(2 * wl / (take_off_distance * 0.85 * density_take_off * 9.80065 * math.pi * oswald_efficiency * aspect_ratio)) + 4 * 11 * 2 / take_off_distance)
        numpy.append(wl, q)
    return T_over_W
    
