import numpy as np

class weight_calculation:
    def __init__(self, n_stringers_up = np.array([]), n_stringers_down = np.array([]), stringer_area=0.0, skin_thickness=0.0, spar_thickness=0.0):
        self.n_stringers_up = n_stringers_up
        self.n_stringers_down = n_stringers_down
        self.stringer_area = stringer_area
        self.skin_thickness = skin_thickness
        self.spar_thickness = spar_thickness
        self.rib_thickness = 123 ## modify
        self.rib_lightening_coeff = 0.123 ## modify, how much is metal
        self.span = 123 ## modify
        self.density = 3000 ## kg/m^3
        self.c_root = 123 ## modify
        self.wing_box_coordinates = [[0.6, 0.0781825],
                                     [0.1, 0.0781825],
                                     [0.1, -0.020718342983016523],
                                     [0.6, -0.026055198062595653]]
        self.taper_ratio = 123 ## modify
        self.front_spar_h_frac = 0
        self.rear_spar_h_frac = 0
        self.top_stringers_c_frac = 0
        self.bottom_stringers_c_frac = 0
        self.airfoil_coordinates = np.array([])

    def read_coordinate_file(self):
        airfoil_coordinates = []
        with open("WP4/NACA64714 a=0.0.dat", "r") as file:
            for line in file:
                ## Skip blank lines or header text
                line = line.strip()
                if not line or not line[0].isdigit():
                    continue

                x, y = map(float, line.split())
                airfoil_coordinates.append([x, y])
        self.airfoil_coordinates = np.array(airfoil_coordinates)

    def chord_length_at_y(self, y_pos):
        semi_span = self.span / 2
        return self.c_root * (1 - (1 - self.taper_ratio) * (y_pos / semi_span))

    def box_fractional_lengths(self):
        self.front_spar_h_frac = self.wing_box_coordinates[1][1] - self.wing_box_coordinates[2][1]
        self.rear_spar_h_frac = self.wing_box_coordinates[0][1] - self.wing_box_coordinates[3][1]
        self.top_stringers_c_frac = self.wing_box_coordinates[0][0] - self.wing_box_coordinates[1][0]
        self.bottom_stringers_c_frac = np.hypot(
            self.wing_box_coordinates[3][0] - self.wing_box_coordinates[2][0],
            self.wing_box_coordinates[3][1] - self.wing_box_coordinates[2][1]
        )
    
    def spar_weight(self, which):
        root_height, tip_height = 0, 0
        if which == "front":
            root_height = self.chord_length_at_y(0) * self.front_spar_h_frac
            tip_height = self.chord_length_at_y(self.span/2) * self.front_spar_h_frac
        elif which == "rear":
            root_height = self.chord_length_at_y(0) * self.rear_spar_h_frac
            tip_height = self.chord_length_at_y(self.span/2) * self.rear_spar_h_frac
        side_area = (root_height + tip_height)/2 * self.span/2
        weight = side_area * self.spar_thickness * self.density
        return weight
    
    def stringers_weight(self, which):
        length = self.span/(2 * 4)
        weights = np.array([])
        if which == "top":
            weights = self.n_stringers_up * length * self.stringer_area * self.density
        elif which == "bottom":
            weights = self.n_stringers_down * length * self.stringer_area * self.density
        return np.sum(weights)

    def airfoil_area_unit_size(self):
        x = self.airfoil_coordinates[:, 0]
        y = self.airfoil_coordinates[:, 1]
        ## Shoelace
        unit_area = 1/2*abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
        return unit_area
    
    def rib_weight(self, y_pos):
        c_at_y = self.chord_length_at_y(y_pos)
        unit_area = self.airfoil_area_unit_size()
        area = c_at_y**2 * unit_area
        weight = self.rib_lightening_coeff * self.rib_thickness * area * self.density
        return weight
    
    def skin_weight(self):
        root_chord_coord_map = self.c_root * self.airfoil_coordinates
        tip_chord_coord_map = self.taper_ratio * self.c_root * self.airfoil_coordinates
        N = root_chord_coord_map.shape[0]
        total_area = 0.0
        for i in range(N-1):
            quad = np.array([
                root_chord_coord_map[i],
                root_chord_coord_map[i+1],
                tip_chord_coord_map[i+1],
                tip_chord_coord_map[i]
            ])
            ## Shoelace
            x = quad[:,0]
            y = quad[:,1]
            area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
            total_area += area
        weight = total_area * self.skin_thickness * self.density
        return weight
    
    def wing_weight(self):
        weight = 0.0
        self.read_coordinate_file()
        self.box_fractional_lengths()
        weight += self.spar_weight("front")
        weight += self.spar_weight("rear")
        weight += self.stringers_weight("top")
        weight += self.stringers_weight("bottom")
        weight += self.rib_weight(self.span/(2 * 2))
        weight += self.skin_weight()
        return weight



design_1 = weight_calculation(np.array([4, 4, 2, 2]),
                              np.array([4, 4, 2, 2]),
                              0.0001,
                              0.002,
                              0.005)

print(design_1.wing_weight())