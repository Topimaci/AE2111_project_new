

m_oe_over_m_MTO_fraction = 0.6077
m_pl_design = 750
m_pl_maxaximum = 1010


def m_MTO_function(m_pl_des, fuel_mass_fraction, m_oe_over_m_MTO_fraction):
    m_MTO = m_pl_des/(1-fuel_mass_fraction-m_oe_over_m_MTO_fraction)
    return m_MTO

def m_oe_function(m_oe_over_m_MTO_fraction, m_MTO):
    m_oe = m_oe_over_m_MTO_fraction * m_MTO
    return m_oe

def m_f_des_function(m_MTO, m_pl_des, m_oe):
    m_f_des = m_MTO - m_pl_des - m_oe
    return m_f_des

def m_f_ferry_function(m_MTO, m_oe):
    m_f_ferry = m_MTO - m_oe
    return m_f_ferry

def m_f_harmonic_function(m_MTO, m_oe, m_pl_max):
    m_f_harmonic = m_MTO - m_oe - m_pl_max
    return m_f_harmonic


