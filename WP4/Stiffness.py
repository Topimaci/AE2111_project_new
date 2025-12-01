import math as m
import sympy as sp
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Integration import x_grid, T_total
from Moment_Diagram import M_vals



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

## print(spar_length(0, 0.101, 2.874, 1.043, 19.585)) ## Spar_location_fraction, y_coordinate, root_chord, tip_chord, span



##______Values________________________________________________________
b = 19.585    # hard coded for now, should probably be pulled from somewhere in the code later on
max_displ = 0.15 * b
max_tip_rotat_deg = 10   # in degrees
max_tip_rotat_rad = m.radians(max_tip_rotat_deg)      # in radians
E = 71 * 10 ** 9    # Young's modulus
G = 27 * 10 ** 9    # Shear modulus

y = sp.symbols("y")
q1, q2, dtheta = sp.symbols('q1 q2 dtheta')


#_______TO BE REPLACED LATER__________________________________________
M_root = 5e6  # N*m
M_y = M_root * (1 - y / (b/2))**2
T_root = 6e5  # N*m, realistic torsion for business jet wingbox
T = T_root * (1 - y / (b/2))**2
y_breaks = np.array([3, 5, 7]) #list of y-positions where the number of stringers decreases, stringer breaks as np.array([...])
stringer_top_num = np.array([5, 4, 3]) #nummber of stringer at the top per interval (that's why it's a list) in np.array([...])
stringer_bottom_num = np.array([5, 4, 3])  #nummber of stringer at the bottom per interval (that's why it's a list) in np.array([...])


#Linear interpolation of the stringers
string_top_interp = interp1d(y_breaks, stringer_top_num, kind="linear",
    fill_value="extrapolate")
string_bottom_interp = interp1d(y_breaks, stringer_bottom_num, kind="linear",
    fill_value="extrapolate")

spar_list = [lambda y: -0.0128 * y + 0.4, 0.3 * b/2, 0.1] #0.5 is how much of the wing span the spar takes, 0.6 is how much of the chord it takes, measured from left side
#spar_list = [lambda y: 0, 0, 0]

def stiffness_distribution(y_pos, h_fs, h_rs, c_upper, c_lower, t, A_string, spar_list, G): #be careful, G is not an input, but still used in this function
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
        i = 0
        while i < len(spar_list) - 1:  # last element may be chord location
            h_spar_func = spar_list[i]
            y_crit = spar_list[i + 1]

            h_spar_y = h_spar_func(y)
            I_spar_y = 1 / 12 * h_spar_y ** 3 * t

            I_step += sp.Piecewise(
                (I_spar_y, y < y_crit),
                (0, True))
            i += 2  # move to next spar
        
        I_total = I_step + I_string_bottom + I_string_top + I_bottom + I_top + I_fs + I_rs
        '''
        if not spar_list[1] == 0:
            alpha = (((1 - spar_list[2]) * c_upper)/spar_list[1]) * y_pos * b/2
            beta = (1 - spar_list[2]) * c_upper - alpha * spar_list[1]
            a = c_upper - (alpha * y_pos * b/2 + beta)
        else:
            a = spar_list[2] * c_upper
        '''
        a = spar_list[2] * c_upper

        w = c_upper - a 
        A_1 = a * w
        A_2 = w * w
        M = sp.Matrix([
        [2*w+2*a, -w, -2*A_1*G*t],
        [-w, 4*w, -2*A_2*G*t],
        [2*A_1, 2*A_2, 0]])

        rhs = sp.Matrix([1, 0, 0])

        solution = M.LUsolve(rhs)

        q1, q2, dtheta_dy = solution
        J = 1 / (G * dtheta_dy)

    else:
        I_total = I_string_bottom + I_string_top + I_bottom + I_top + I_fs + I_rs
        A = (h_fs + h_rs) / 2 * c_upper
        circ = 1/ t * (h_fs + c_upper + h_rs + c_lower)
        J = 4 * A ** 2 / circ
    return I_total, J


#______THIS IS WHERE WE CALL THE FUNCTION, ALL OF THE VALUES MUST BE REPLACED WITH THE CORRECT ONES
I_xx, J = stiffness_distribution(x_grid, 0.3, 1, 1.2, 0.05, 0.2, 0.1, spar_list, G)

d2v_dy2 = -M_vals / (E * I_xx)
dth_dy  =  T_total   / (G * J)

f_d2v = sp.lambdify(y, d2v_dy2, "numpy")
f_th  = sp.lambdify(y, dth_dy,   "numpy")

# First integration
estimate_dv, error_dv = integrate.quad(f_d2v, 0, b)

# Nested integration for displacement
def slope(Y):
    return integrate.quad(f_d2v, 0, Y)[0]

# Deflection
estimate_v, error_v = integrate.quad(slope, 0, b/2)

# Twist
estimate_th, error_th = integrate.quad(f_th, 0, b/2)

'''
print(estimate_v, error_v)
print(estimate_th, error_th)
print(m.degrees(estimate_th))
print("I_xx:", I_xx)
print("J:", J)
'''


##______Output results________________________________________________________
## Get the box coordinates for spar length calculations
def diagram_plotter(spar_location_fraction1, spar_location_fraction2, root_chord, tip_chord, b, t, A_string, spar_list, steps):
    spar1_coor1, spar1_coor2, spar1_coor3, spar1_coor4 = spar_position(Airfoil_coordinates, spar_location_fraction1)
    spar2_coor1, spar2_coor2, spar2_coor3, spar2_coor4 = spar_position(Airfoil_coordinates, spar_location_fraction2)
    y_at_top_stringer = top_stringer_y_coord(spar1_coor1, spar1_coor2, spar2_coor1, spar2_coor2, spar_location_fraction1)
    y1_at_bot_stringer, y2_at_bot_stringer = bot_stringer_y_coords(spar1_coor3, spar1_coor4, spar2_coor3, spar2_coor4, spar_location_fraction1, spar_location_fraction2)
    box_coordinates = [[spar_location_fraction2, y_at_top_stringer], [spar_location_fraction1, y_at_top_stringer],
                       [spar_location_fraction1, y1_at_bot_stringer], [spar_location_fraction2, y2_at_bot_stringer]]

    y_array = np.linspace(0, b/2, steps)
    dv_dy2_array = np.zeros_like(y_array)
    dtheta_dy_array = np.zeros_like(y_array)

    # Compute d2v/dy2 and dtheta/dy numerically
    for i, y_i in enumerate(y_array):
        front_spar_length = spar_length(spar_location_fraction1, y_i, root_chord, tip_chord, b, spar_location_fraction1, spar_location_fraction2)
        rear_spar_length = spar_length(spar_location_fraction2, y_i, root_chord, tip_chord, b, spar_location_fraction1, spar_location_fraction2)
        h_fs, h_rs = front_spar_length, rear_spar_length

        a = (tip_chord - root_chord)/(b/2)
        chord_length_at_y = root_chord + a * y_i
        c_upper = abs(box_coordinates[0][0] - box_coordinates[1][0]) * chord_length_at_y
        c_lower = m.sqrt((box_coordinates[2][0] - box_coordinates[3][0])**2 + (box_coordinates[2][1] - box_coordinates[3][1])**2) * chord_length_at_y

        I_xx_sym, J = stiffness_distribution(y_i, h_fs, h_rs, c_upper, c_lower, t, A_string, spar_list)
        I_xx = float(I_xx_sym.subs(y, y_i))

        # Make numeric functions from symbolic expressions
        f_M = sp.lambdify(y, M_y, "numpy")
        f_T = sp.lambdify(y, T, "numpy")

        # Inside the loop over y_array
        dv_dy2_array[i]   = - f_M(y_i) / (E * I_xx)
        dtheta_dy_array[i] = f_T(y_i) / (G * J)
    

    # Numerical integration using cumulative trapezoid
    v_array = np.cumsum((dv_dy2_array[:-1] + dv_dy2_array[1:]) / 2 * np.diff(y_array))
    theta_array = np.cumsum((dtheta_dy_array[:-1] + dtheta_dy_array[1:]) / 2 * np.diff(y_array))

    # Add zero at the root
    v_array = np.insert(v_array, 0, 0)
    theta_array = np.insert(theta_array, 0, 0)

    plt.plot(y_array, v_array, label="Deflection (m)")
    plt.xlabel("Spanwise position (m)")
    plt.ylabel("Deflection")
    plt.show()
    plt.plot(y_array, np.degrees(theta_array), label="Twist (deg)")
    plt.xlabel("Spanwise position (m)")
    plt.ylabel("Twist")
    plt.vlines(b/2, min(np.degrees(theta_array)), max(np.degrees(theta_array)), colors='r', linestyles='dashed', label='Wing Tip')
    plt.show()


diagram_plotter(0.1, 0.6, 2.874, 1.043, b, 0.003, 0.2, spar_list, 60)