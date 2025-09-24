import numpy as np
import math as m

#__wing loadings, mass fraction, climb gradient, cd_0, oswald number for current config, AR for current config, density, cl_max and bypass ratio is the input

def climb_grad(wing_loadings, mass_fraction, climb_gradient, cd_0, oswald_number, aspect_ratio, density, cl_max, B):
    
    t_over_w = np.array([])
    for wing_loading in wing_loadings:
        velocity = np.sqrt(2*wing_loading/(density*cl_max))
        mach = velocity/np.sqrt(287*1.4*288.15)
        total_temps = 288.15*(1+0.4/2*mach**2)/288.15
        total_pressure = 101325*(1+0.4/2*mach**2)**3.5/101325
        laps_rate = total_pressure*(1-(0.43+0.014*B)*m.sqrt(mach))
        tw = mass_fraction/laps_rate * (climb_gradient/100+(2*m.sqrt(cd_0/(m.pi * aspect_ratio * oswald_number))))
        t_over_w = np.append(t_over_w,tw)
    

    return t_over_w

    
