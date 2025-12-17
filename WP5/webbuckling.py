#Packages
import math as m
import numpy as np

#Front spar and Rear spar files

FRONT_SPAR_FILE = "h_front_spar.npy"
REAR_SPAR_FILE = "h_rear_spar.npy"
h_fs = np.load(FRONT_SPAR_FILE)
h_rs = np.load(REAR_SPAR_FILE)
print(h_fs[0])
print(h_rs[0])

#Material properties
Pois = 0.33
E = 72.4 * 10**9 #Pa
b = 1 #m
t = 0.06 #m

k_v = 1.3

def average_shear_stress(V, h_s, h_r, t_f, t_r):
    return V / (h_s * t_f + h_r * t_r)

def critical_shear_stress(k_s):
    return m.pi ** 2 * k_s * E / (12 - 12 * Pois ** 2) * (t/b)

def max_shear_stress(ave_shear_stress):
    return k_v * ave_shear_stress
