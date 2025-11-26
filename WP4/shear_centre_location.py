import math
from scipy.integrate import quad

#Assumptions:
#shear modulus G is constant through the section
#horizontal axis is symmetry axis

chord   = 1000
x_front = 0.33 * chord
x_rear  = 0.75 * chord
c_front = 0.1205 * chord
c_rear  = 0.0783 * chord

t_front = 10
t_rear  = 10
t_web   = 10


#input to check if code is correct
#chord = 1
# x_front = 0
# x_rear  = 786
# c_front = 600
# c_rear  = 300
#
# t_front = 12
# t_rear  = 8
# t_web   = 10


# Diagonal web length
a = math.sqrt((x_rear - x_front)**2 + ((c_front - c_rear)/2)**2)

#Second moment of area of the section about x-axis is given by
def inertia_moment_xx():
    def web_integrand(s):
        y = (c_rear/2) + (c_rear/(2*a)) * s
        return t_web * y**2

    integral_web = quad(web_integrand, 0, a)
    Ixx = (1/12) * t_front * c_front**3 + (1/12) * t_rear * c_rear**3 + 2 * integral_web[0]
    return Ixx

#We now obtain the q_b shear flow distribution by 'cutting' the beam section at the mid-point O of the wall CB. Thus, since y=s_a we have
def shear_flow_ob(v_y, Ixx, s_a):
    return -v_y / Ixx * t_rear * s_a**2 / 2

#For the wall BA where y = c_rear/2 + c_rear * s_b / (2a)
def shear_flow_ba(v_y, Ixx, s_b):
    return -v_y / Ixx * (t_web * c_rear * (s_b/2 + s_b**2 / (4*a)) + 0.125 * t_rear * c_rear**2)

#In the wall AD, y = c_front/2 - s_c so that
def shear_flow_ad(v_y, Ixx, s_c):
    return -v_y / Ixx * (t_front * (c_front * s_c /2 - s_c**2 / 2) + (3/4 * t_web * c_rear * a + 0.125 * t_rear * c_rear**2))


def shear_flow_const(v_y):
    Ixx = inertia_moment_xx()

    line_integral = c_front / t_front + 2 * a / t_web + c_rear / t_rear

    #integral along OB
    def integrand_ob(s_a, v_y, Ixx):
        q_b = shear_flow_ob(v_y, Ixx, s_a)
        return q_b / t_rear

    integral_ob = quad(integrand_ob, 0, c_rear / 2, args=(v_y, Ixx))[0]

    #integral along BA
    def integrand_ba(s_b, v_y, Ixx):
        q_b = shear_flow_ba(v_y, Ixx, s_b)
        return q_b / t_web

    integral_ba = quad(integrand_ba, 0, a, args=(v_y, Ixx))[0]

    #integral along AD
    def integrand_ad(s_c, v_y, Ixx):
        q_b = shear_flow_ad(v_y, Ixx, s_c)
        return q_b / t_front

    integral_ad = quad(integrand_ad, 0, c_front / 2, args=(v_y, Ixx))[0]

    return -2 / line_integral * (integral_ad + integral_ba + integral_ob)

def integrate_shear_flow_ob(v_y, Ixx, q_s):
    def integrand(s_a, v_y, Ixx, q_s):
        moment_arm = x_rear - x_front
        return (shear_flow_ob(v_y, Ixx, s_a) + q_s) * moment_arm
    result = quad(integrand, 0, c_rear/2, args=(v_y, Ixx, q_s))
    return result[0]

def integrate_shear_flow_ba(v_y, Ixx, q_s):
    def integrand(s_b, v_y, Ixx, q_s):
        moment_arm = c_front/2 * (x_rear - x_front)/a
        return (shear_flow_ba(v_y, Ixx, s_b) + q_s) * moment_arm
    result = quad(integrand, 0, a, args=(v_y, Ixx, q_s))
    return result[0]

# def integrate_shear_flow_ad(v_y, Ixx, q_s):
#     def integrand(s_c, v_y, Ixx, q_s):
#         return (shear_flow_ad(v_y, Ixx, s_c) + q_s) * (x_rear - x_front)
#     result,_ = quad(integrand, 0, c_front/2, args=(v_y, Ixx, q_s))
#     return result


def shear_center_non_dim(v_y=1.0):
    Ixx = inertia_moment_xx()
    q_s0 = shear_flow_const(v_y)

    M_ob = integrate_shear_flow_ob(v_y, Ixx, q_s0)
    M_ba = integrate_shear_flow_ba(v_y, Ixx, q_s0)

    ksi = 2 * (M_ob + M_ba) / v_y
    return (ksi + x_front) / chord
