import WP2.master_range_mass as mrm
import variables.dynamic_variables as dv

mass_to_new = mrm.m_MTO
mass_oe_new = mrm.m_oe
mass_fuel_new = mrm.m_f_des

S_wing = dv.S_w
S_wing_new = S_wing

while S_wing_new/S_wing <= 0.05:

    mass_to = mass_to_new
    mass_oe = mass_oe_new
    mass_fuel = mass_fuel_new

    
