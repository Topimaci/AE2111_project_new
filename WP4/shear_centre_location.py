import math
from scipy.integrate import quad

chord   = 2.9
x_front = 0.2 * chord
x_rear  = 0.7 * chord
c_front = 0.121 * chord
c_rear  = 0.0787 * chord

t_front = 0.01
t_rear  = 0.01
t_web   = 0.01

# Diagonal web length
a = math.sqrt((x_rear - x_front)**2 + ((c_front - c_rear)/2)**2)

def inertia_moment_xx():
    def web_integrand(s):
        y = (c_rear/2) + (c_rear/(2*a)) * s
        return t_web * y**2

    integral_web = quad(web_integrand, 0, a)
    Ixx = (1/12) * t_front * c_front**3 + (1/12) * t_rear * c_rear**3 + 2 * integral_web[0]
    return Ixx


def shear_flow_ob(v_y, Ixx, s_a):
    return -v_y / Ixx * t_rear * s_a**2 / 2

def shear_flow_ba(v_y, Ixx, s_b):
    return -v_y / Ixx * (t_web * c_rear * (s_b/2 + s_b**2 / (4*a)) + 0.25 * t_rear * c_rear**2)

def shear_flow_ad(v_y, Ixx, s_c):
    return -v_y / Ixx * (t_front * (c_front * s_c /2 - s_c**2 / 2) + (3/4 * t_web * c_rear * a + 0.25 * t_rear * c_rear**2))


def shear_flow_const(v_y):
    Ixx = inertia_moment_xx()
    line_integral = c_front / t_front + 2 * a / t_web + c_rear / t_rear

    q_s0 = -2 / line_integral * (
        shear_flow_ob(v_y, Ixx, c_rear/2) / t_rear +
        shear_flow_ba(v_y, Ixx, a) / t_web +
        shear_flow_ad(v_y, Ixx, c_front/2) / t_front
    )
    return q_s0


def integrate_shear_flow_ob(v_y, Ixx, q_s):
    def integrand(s_a, v_y, Ixx, q_s):
        return (shear_flow_ob(v_y, Ixx, s_a) + q_s) * (x_rear - x_front)
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

    ksi = (M_ob + M_ba) / v_y
    return ksi / chord  # non-dimensional


print("Non-dimensional shear-center position ksi =", shear_center_non_dim(3))