import math

def fuel_volume(root_chord, taper, span, m_fuel):
    beginning_fuel_tank = 0
    end_fuel_tank = 0.9
    t_over_c = 0.14
    LE_pos = 0.18
    TE_pos = 0.58

    Landing_gear_beginning = 0.3
    Landing_gear_end = 0.4

    distance_tank1 = span/2 *beginning_fuel_tank - span/2 *Landing_gear_beginning
    distance_tank2 = span/2 *Landing_gear_end - span/2 *end_fuel_tank

    chord_beginning = root_chord - root_chord*(1-taper)*beginning_fuel_tank
    chord_landing_gear_beginning = root_chord - root_chord*(1-taper)*Landing_gear_beginning
    chord_landing_gear_end = root_chord - root_chord*(1-taper)*Landing_gear_end
    chord_end = root_chord - root_chord*(1-taper)*end_fuel_tank


    Thickness_LE_beginning = 2*chord_beginning *(t_over_c/2)*(0.2969*math.sqrt(LE_pos)-0.126*LE_pos-0.3516*LE_pos**2 + 0.2843*LE_pos**3 - 0.1015*LE_pos**4)
    Thickness_TE_beginning = 2*chord_beginning *(t_over_c/2)*(0.2969*math.sqrt(TE_pos)-0.126*TE_pos-0.3516*TE_pos**2 + 0.2843*TE_pos**3 - 0.1015*TE_pos**4)
    average_thickness_beginning = (Thickness_LE_beginning+Thickness_TE_beginning)/2
    Area_beginning = (TE_pos - LE_pos)*chord_beginning*average_thickness_beginning

    Thickness_LE_LG_beginning = 2*chord_landing_gear_beginning *(t_over_c/2)*(0.2969*math.sqrt(LE_pos)-0.126*LE_pos-0.3516*LE_pos**2 + 0.2843*LE_pos**3 - 0.1015*LE_pos**4)
    Thickness_TE_LG_beginning = 2*chord_landing_gear_beginning *(t_over_c/2)*(0.2969*math.sqrt(TE_pos)-0.126*TE_pos-0.3516*TE_pos**2 + 0.2843*TE_pos**3 - 0.1015*TE_pos**4)
    average_thickness_LG_beginning = (Thickness_LE_LG_beginning+Thickness_TE_LG_beginning)/2
    Area_LG_beginning = (TE_pos - LE_pos)*chord_landing_gear_beginning*average_thickness_LG_beginning

    Volume_tank1 = distance_tank1/3 *(Area_beginning+Area_LG_beginning +math.sqrt(Area_beginning*Area_LG_beginning))
    

    Thickness_LE_LG_end = 2*chord_landing_gear_end *(t_over_c/2)*(0.2969*math.sqrt(LE_pos)-0.126*LE_pos-0.3516*LE_pos**2 + 0.2843*LE_pos**3 - 0.1015*LE_pos**4)
    Thickness_TE_LG_end = 2*chord_landing_gear_end *(t_over_c/2)*(0.2969*math.sqrt(TE_pos)-0.126*TE_pos-0.3516*TE_pos**2 + 0.2843*TE_pos**3 - 0.1015*TE_pos**4)
    average_thickness_LG_end = (Thickness_LE_LG_end+Thickness_TE_LG_end)/2
    Area_LG_end= (TE_pos - LE_pos)*chord_landing_gear_end*average_thickness_LG_end

    Thickness_LE_end = 2*chord_end *(t_over_c/2)*(0.2969*math.sqrt(LE_pos)-0.126*LE_pos-0.3516*LE_pos**2 + 0.2843*LE_pos**3 - 0.1015*LE_pos**4)
    Thickness_TE_end = 2*chord_end *(t_over_c/2)*(0.2969*math.sqrt(TE_pos)-0.126*TE_pos-0.3516*TE_pos**2 + 0.2843*TE_pos**3 - 0.1015*TE_pos**4)
    average_thickness_end = (Thickness_LE_end+Thickness_TE_end)/2
    Area_end = (TE_pos - LE_pos)*chord_end*average_thickness_end

    Volume_tank2 = distance_tank2/3 *(Area_LG_end + Area_end + math.sqrt(Area_LG_end*Area_end))

    Total_volume_wing = (Volume_tank1 + Volume_tank2) * 2
    density = 800   
    Total_volume_fuel_needed = m_fuel/density

    percentage_fuel_in_wing = Total_volume_wing/Total_volume_fuel_needed
    fuel_mass_in_wing = percentage_fuel_in_wing*m_fuel
    return Total_volume_wing, percentage_fuel_in_wing, m_fuel












