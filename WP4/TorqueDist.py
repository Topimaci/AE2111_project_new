import numpy as np
import scipy as sp

#Torque distribution module for WP4, first integral of the distributed torque function

#Testing scipy integrationmethod

def f(x):
    return np.sqrt(1-(x/40)**2)

estimatef, errorf = sp.integrate.quad(f, 0, 40)

def g(x):
    res = 20 if x < 5 else 0
    return res

estimateg, errorg = sp.integrate.quad(g, 0, 10)

print("Integral of f from 0 to 40 is:", estimatef, "with error estimate:", errorf)
print("Integral of g from 0 to 10 is:", estimateg, "with error estimate:", errorg)