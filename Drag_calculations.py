import math

S_wet_over_S_w = 6
C_f = 0.001
psi = 
phi = 
AR = 



def c_d0_function(S_wet_over_S_w, C_f):
    c_d0 = (C_f * S_wet_over_S_w)
    return c_d0



def e_function(psi, phi, AR):
    e = 1/(math.pi*AR*psi + phi)
    return e


def L_over_D_max_function(AR,e,C_D_0):  #### watch out what kind of C_L you are using!!
    C_L_for_max_L_over_D = math.sqrt(math.pi*AR*e*C_D_0)
    C_D_for_max_L_over_D = 2*C_D_0
    L_over_D_max = C_L_for_max_L_over_D/C_D_for_max_L_over_D

    return  L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D

def c_d_function(C_D_0, C_L, AR, e):
    C_D = C_D_0 + C_L**2/(math.pi*AR*e)
    return C_D


C_D_0_initial = c_d0_function(S_wet_over_S_w, C_f)
e_initial = e_function(psi, phi, AR)

L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D = L_over_D_max_function(AR,e_initial,C_D_0_initial)