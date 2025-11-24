import math as m
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

##______Spar length based on airfoil and y-position________________________________________________________
Airfoil_coordinates = []

## Read airfoil coordinates from file
with open("WP4/NACA64714 a=0.0.dat", "r") as file:
    for line in file:
        ## Skip blank lines or header text
        line = line.strip()
        if not line or not line[0].isdigit():
            continue

        x, y = map(float, line.split())
        Airfoil_coordinates.append([x, y])


## For a single spar calculate the top and bottom closest points from the airfoil map
def spar_position(Airfoil_coordinates, spar_location_fraction):
    closest_coordinates_1 = [0, 0]
    closest_coordinates_2 = [0, 0]
    closest_coordinates_3 = [0, 0]
    closest_coordinates_4 = [0, 0]
    previous_coords = [0, 0]
    one_found = True

    ## Find the boundary coordinates for the spar location
    for coords in Airfoil_coordinates:
        ## Top spar coordinates
        if spar_location_fraction > coords[0] and one_found:
            closest_coordinates_2 = coords
            closest_coordinates_1 = previous_coords
            one_found = False
        ## Bottom spar coordinates
        if spar_location_fraction < coords[0] and not one_found:
            closest_coordinates_4 = coords
            closest_coordinates_3 = previous_coords
            break
        previous_coords = coords
    return(closest_coordinates_1, closest_coordinates_2, closest_coordinates_3, closest_coordinates_4)

## Find the top stringer y coordinates based on the front and rear spar, they are horizontal
def top_stringer_y_coord(spar1_coor1, spar1_coor2, spar2_coor1, spar2_coor2, spar_location_fraction1):
    ## Find the bottom of the top stringer
    y_top_spar1 = max(spar1_coor1[1], spar1_coor2[1])
    y_top_spar2 = max(spar2_coor1[1], spar2_coor2[1])
    top_y = min(y_top_spar1, y_top_spar2)
    ## Find the top of the bottom stringer
    y_bot_spar1 = min(spar1_coor1[1], spar1_coor2[1])
    y_bot_spar2 = min(spar2_coor1[1], spar2_coor2[1])
    bot_y = max(y_bot_spar1, y_bot_spar2)
    
    if top_y - bot_y <= 0: ## If there is no overlap, take the average
        return 1/2 * (y_top_spar1 + y_top_spar2)
    else: ## If there is overlap, find the y coordinate at the spar chord-wise location via a line between two of the bounding points
        a = (spar1_coor2[1] - spar1_coor1[1])/(spar1_coor2[0] - spar1_coor1[0])
        b = spar1_coor1[1] - a * spar1_coor1[0]
        y_at_spar_location = a * spar_location_fraction1 + b
        return y_at_spar_location

## Find the bottom stringer y coordinates based on the front and rear spar
def bot_stringer_y_coords(spar1_coor3, spar1_coor4, spar2_coor3, spar2_coor4, spar_location_fraction1, spar_location_fraction2):
    a1 = (spar1_coor4[1] - spar1_coor3[1])/(spar1_coor4[0] - spar1_coor3[0])
    b1 = spar1_coor3[1] - a1 * spar1_coor3[0]
    y1_at_spar_location = a1 * spar_location_fraction1 + b1

    a2 = (spar2_coor4[1] - spar2_coor3[1])/(spar2_coor4[0] - spar2_coor3[0])
    b2 = spar2_coor3[1] - a2 * spar2_coor3[0]
    y2_at_spar_location = a2 * spar_location_fraction2 + b2

    return y1_at_spar_location, y2_at_spar_location

## Find the spar length fraction based on the box coordinates
def spar_length_fraction(box_coordinates, spar_location_fraction):
    a = (box_coordinates[3][1] - box_coordinates[2][1])/(box_coordinates[3][0] - box_coordinates[2][0])
    b = box_coordinates[2][1] - a * box_coordinates[2][0]
    y_at_spar_location = a * spar_location_fraction + b
    fraction = box_coordinates[0][1] - y_at_spar_location
    return fraction

## Final spar length calculation function
def spar_length(spar_location_fraction, y_coordinate, root_chord, tip_chord, span, spar_location_fraction1, spar_location_fraction2):
    spar1_coor1, spar1_coor2, spar1_coor3, spar1_coor4 = spar_position(Airfoil_coordinates, spar_location_fraction1) ## Front spar, top right, top left, bottom left, bottom right
    spar2_coor1, spar2_coor2, spar2_coor3, spar2_coor4 = spar_position(Airfoil_coordinates, spar_location_fraction2) ## Rear spar, top right, top left, bottom left, bottom right
    y_at_top_stringer = top_stringer_y_coord(spar1_coor1, spar1_coor2, spar2_coor1, spar2_coor2, spar_location_fraction1)
    y1_at_bot_stringer, y2_at_bot_stringer = bot_stringer_y_coords(spar1_coor3, spar1_coor4, spar2_coor3, spar2_coor4, spar_location_fraction1, spar_location_fraction2)
    box_coordinates = [[spar_location_fraction2, y_at_top_stringer], [spar_location_fraction1, y_at_top_stringer],
                    [spar_location_fraction1, y1_at_bot_stringer], [spar_location_fraction2, y2_at_bot_stringer]] ## [top right, top left, bottom left, bottom right]
    length_fraction = spar_length_fraction(box_coordinates, spar_location_fraction)
    a = (tip_chord - root_chord)/(span/2)
    chord_length_at_y = root_chord + a * y_coordinate
    spar_length = length_fraction * chord_length_at_y

    ## plot the airfoil and torsion box
    '''
    xpoints = []
    ypoints = []
    for coords in Airfoil_coordinates:
        xpoints.append(coords[0])
        ypoints.append(coords[1])

    xbox = []
    ybox = []
    for coords in box_coordinates:
        xbox.append(coords[0])
        ybox.append(coords[1])

    plt.axis('equal')
    plt.plot(xpoints, ypoints)
    #plt.plot([spar1_coor1[0], spar1_coor2[0], spar1_coor3[0], spar1_coor4[0]],
    #        [spar1_coor1[1], spar1_coor2[1], spar1_coor3[1], spar1_coor4[1]], 'o', color='red')
    #plt.plot([spar2_coor1[0], spar2_coor2[0], spar2_coor3[0], spar2_coor4[0]],
    #        [spar2_coor1[1], spar2_coor2[1], spar2_coor3[1], spar2_coor4[1]], 'o', color='red')
    plt.plot(xbox + [xbox[0]], ybox + [ybox[0]], color='green')
    plt.show()
    '''

    return spar_length

print(spar_length(0, 0.101, 2.874, 1.043, 19.585)) ## Spar_location_fraction, y_coordinate, root_chord, tip_chord, span



##______Deflection & twist calculation________________________________________________________
b = 19.585    # hard coded for now, should probably be pulled from somewhere in the code later on
max_displ = 0.15 * b
max_tip_rotat_deg = 10   # in degrees
max_tip_rotat_rad = m.radians(max_tip_rotat_deg)      # in radians
E = 71 * 10 ** 9    # Young's modulus
G = 27 * 10 ** 9    # Shear modulus

M_x =   # Import Moment function of M(y)
T =     # Import torque distribution function

y_breaks = #stringer breaks as np.array([...])
stringer_top_num = #nummber of stringer on these intervals np.array([...])
stringer_bottom_num = #same thing

#Linear interpolation of the stringers
string_top_interp = interp1d(y_breaks, stringer_top_num, kind="linear",
    fill_value="extrapolate")
string_bottom_interp = interp1d(y_breaks, stringer_bottom_num, kind="linear",
    fill_value="extrapolate")

y = sp.symbols("y")

spar_list = [lambda y: 0.4 * y + 0.1, 0.5, 1] # functions should be replaced, this is just an example, 0.5 is how much of the wing span the spar takes, 1 is how much of the chord it takes, measured from left side


def stiffness_distribution(y_pos, h_fs, h_rs, c_upper, c_lower, t, A_string, spar_list):
    # I Moment of Inertia Calculations
    #neutral axis
    x_c = (h_rs ** 2 + h_fs ** 2 + h_fs * h_rs) / (3 * (h_rs + h_fs))
    #Spar inertias
    I_fs = 1/12 * h_fs ** 3 * t + h_fs * t * (x_c - h_fs/2)**2
    I_rs = 1/12 * h_rs ** 3 * t + h_rs * t * (x_c - h_rs/2)**2
    #Skin inertias
    I_top = t * c_upper * (t/2 - x_c)**2
    I_bottom = t * c_lower * (((h_fs - x_c) + (h_rs - x_c))/2) ** 2

    num_top = float(string_top_interp(y_pos))
    num_bottom = float(string_bottom_interp(y_pos))

    #stringer inertias
    I_string_top = (A_string * (t - x_c)**2) * num_top
    I_string_bottom = (A_string * (((h_fs - x_c) + (h_rs - x_c))/2) ** 2) * num_bottom
       
    if spar_list != []:
        I_step = 0
        for h_spar_func, y_crit in spar_list:

            # h_spar is a sympy expression of y
            h_spar_y = h_spar_func(y)
            # compute I_spar(y)
            I_spar_y = 1/12 * h_spar_y**3 * t
            # add step contribution
            I_step += sp.Piecewise(
                (I_spar_y, y < y_crit),
                (0, True)
            )
        I_total = I_step + I_string_bottom + I_string_top + I_bottom + I_top + I_fs + I_rs
        a = spar_list[2]
        w = c_upper - a 
        lefthand_matrix = np.array([[(2*w+2*a), -w, -2*a*w*G*t], [-w, 4*b, - 2*w**2*G*t], [2*a*w, 2 * w**2, 0]])
        righthand_matrix = np.array([0, 0, 1])
        solution = np.linalg.solve(lefthand_matrix, righthand_matrix)
        q1, q2, dtheta_dy = solution
        J = 1 / (G * dtheta_dy)

    else:
        I_total = I_string_bottom + I_string_top + I_bottom + I_top + I_fs + I_rs
        A = (h_fs + h_rs) / 2 * c_upper
        circ = 1/ t * (h_fs + c_upper + h_rs + c_lower)
        J = 4 * A ** 2 / circ
    return I_total, J

I_xx, J = stiffness_distribution()
d2v_dy2 = - M_x / (E * I_xx)
dth_dy = T / (G * J)

estimate_dv, error_dv = sp.integrate.quad(d2v_dy2, 0, b)
estimate_v, error_v = sp.integrate.quad(estimate_dv, 0, b)
estimate_th, error_th = sp.integrate.quad(dth_dy, 0, b)


##______Output results________________________________________________________
## Get the box coordinates for spar length calculations
def diagram_plotter(spar_location_fraction1, spar_location_fraction2, root_chord, tip_chord, b, t, A_string, spar_list):
    spar1_coor1, spar1_coor2, spar1_coor3, spar1_coor4 = spar_position(Airfoil_coordinates, spar_location_fraction1) ## Front spar, top right, top left, bottom left, bottom right
    spar2_coor1, spar2_coor2, spar2_coor3, spar2_coor4 = spar_position(Airfoil_coordinates, spar_location_fraction2) ## Rear spar, top right, top left, bottom left, bottom right
    y_at_top_stringer = top_stringer_y_coord(spar1_coor1, spar1_coor2, spar2_coor1, spar2_coor2, spar_location_fraction1)
    y1_at_bot_stringer, y2_at_bot_stringer = bot_stringer_y_coords(spar1_coor3, spar1_coor4, spar2_coor3, spar2_coor4, spar_location_fraction1, spar_location_fraction2)
    box_coordinates = [[spar_location_fraction2, y_at_top_stringer], [spar_location_fraction1, y_at_top_stringer],
                    [spar_location_fraction1, y1_at_bot_stringer], [spar_location_fraction2, y2_at_bot_stringer]] ## [top right, top left, bottom left, bottom right]

    deflectionY = []
    deflectionZ = []
    twistY = []
    twist_deg = []

    for i in range(0, 1/(b/2), 0.01):
        front_spar_length = spar_length(spar_location_fraction1, i * (b/2), 2.874, 1.043, b)
        rear_spar_length = spar_length(spar_location_fraction2, i * (b/2), 2.874, 1.043, b)
        h_fs = front_spar_length
        h_rs = rear_spar_length

        a = (tip_chord - root_chord)/(b/2)
        chord_length_at_y = root_chord + a * i * (b/2)
        c_upper = abs(box_coordinates[0][0] - box_coordinates[1][0]) * chord_length_at_y
        c_lower = m.sqrt((box_coordinates[2][0] - box_coordinates[3][0])**2 + (box_coordinates[2][1] - box_coordinates[3][1])**2) * chord_length_at_y

        I_xx, J = stiffness_distribution(i, h_fs, h_rs, c_upper, c_lower, t, A_string, spar_list)
        d2v_dy2 = - M_x / (E * I_xx)
        dth_dy = T / (G * J)

        estimate_dv, error_dv = sp.integrate.quad(d2v_dy2, 0, b)
        estimate_v, error_v = sp.integrate.quad(estimate_dv, 0, b)
        estimate_th, error_th = sp.integrate.quad(dth_dy, 0, b)

        deflectionY.append(i*b/2)
        twistY.append(i*b/2)
        deflectionZ.append(estimate_v)
        twist_deg.append(m.degrees(estimate_th))
    
    plt.plot(deflectionY, deflectionZ)
    plt.xlabel("Spanwise position (m)")
    plt.ylabel("Deflection (m)")
    plt.show()


diagram_plotter(0.1, 0.6, 2.874, 1.043, b, 0.003, 0.0001, spar_list)