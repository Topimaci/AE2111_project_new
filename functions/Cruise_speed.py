import numpy
import math

#### thrust lapse = 0.24  beta = 0.95


def cruise_speed(beta, thrust_lapse, Wing_loading, C_D_0, density, v_cr, AR, e):
    T_over_W = numpy.array([])
    for i in Wing_loading:
        x = (beta/thrust_lapse)*((C_D_0*0.5*density*v_cr**2)/(beta*i)+(beta*i)/(math.pi*AR*e*0.5*density*v_cr**2))
        T_over_W = numpy.append(T_over_W, x)

    return T_over_W

