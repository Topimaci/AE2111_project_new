import math as m

m_MTO = 12520 ## kg
wing_span = 17.45
roll_rate = m.radians(20/5) ## rad/s
density_ISO = 1.225225
V_stall = 56 ## m/s
V_MC = V_stall * 1.3
C_lda = 0.8
max_aileron_deflection = m.radians(20)


def aileron_geometry_calculation():
    I_x_longitudinal = 1/4 * m_MTO * wing_span**2
    rolling_moment = I_x_longitudinal * roll_rate
    dynamic_pressure = 0.5 * density_ISO * (V_MC)**2
    S_a = rolling_moment / (C_lda * dynamic_pressure * wing_span * max_aileron_deflection)

    return S_a

print(aileron_geometry_calculation())
