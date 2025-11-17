import variables.dynamic_variables as dv
import functions.Class_II_weight_estimations as c2w
import WP2.main as ma
import WP3.main as wp3
import functions.Planform_DESIGN1 as pd
import numpy as np
import functions.Drag_calculations_class_II as D2
import functions.Empennage_sizing as es
import functions.Engine_types as et
import WP2.HLDs as hld
import functions.Minimum_speed as ms
import functions.Climb_rate as cr
import functions.Cruise_speed as cs
import functions.Landing_field_length as lfl
import functions.climb_grad as cg
import functions.Take_off_distance as td
import variables.fixed_values as fv
import functions.Fuel_Volume as fuelv
import functions.Range_calculations as rc
import pandas as pa
import matplotlib.pyplot as plt
import math as math



rho = 1.225 
W = 14266 
n = 1.5 #Loading factor 

#Speeds 
V_s0 = math.sqrt(2*W/(rho*wp3.S_wing*fv.C_L_max_landing)) #Stall speed with flaps extended
V_s1 = math.sqrt(2*W/(rho*wp3.S_wing*fv.C_L_max_landing)) #Stall speed with flaps retracted 
V_A = V_s1 * math.sqrt(n) #Manoeuvring speed 
V_C = fv.v_cr #Design cruise speed 
#Design dive speed 


 #Design wing-flap speed 

n_max = 2.1 + 24000/(W + 10000)


V = int()
n = (V/V_s1)^2

plt.plot(V,n)
plt.show()

