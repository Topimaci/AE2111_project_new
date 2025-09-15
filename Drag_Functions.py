import math

def C_D_function(C_D_0, C_L, AR, e):
    C_D = C_D_0 + C_L**2/(math.pi*AR*e)
    return C_D


def L_over_D_max_function(AR,e,C_D_0):  #### watch out what kind of C_L you are using!!
    C_L_for_max_L_over_D = math.sqrt(math.pi*AR*e*C_D_0)
    C_D_for_max_L_over_D = 2*C_D_0
    L_over_D_max = C_L_for_max_L_over_D/C_D_for_max_L_over_D

    return  L_over_D_max, C_L_for_max_L_over_D, C_D_for_max_L_over_D

def C_D_0_function(C_f_over_S_wet, S_w):
    C_D_0 = (C_f * C_f_over_S_wet)
    return C_D_0

def e_function(psi, phi, AR):
    e = 1/(math.pi*AR*psi + phi)
    return e

