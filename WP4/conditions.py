#load_factor = 3 #+LC4
#load_factor = -1.5 #-LC4 
#load_factor = -4 #+LC9
#load_factor = -1.5 #-LC9
#load_factor = 4 #+LC15
#load_factor = -1.5  #-LC15
#load_factor = 3 #+LC22
#load_factor = -1.5 #-LC22
#load_factor = 3 #+LC28
#load_factor = -1.5 #-LC28
#velocity = 60.31   #+LC4
#velocity = 107.91 #+LC9
#velocity = 86.52 #+LC15
#velocity = 86.3 #+LC22
#velocity = 69.19 #+LC28
#weight = 140000 #+LC4
#weight = 140000 #+LC9
#weight = 84669 #+LC15
#weight = 140000 #+LC22
#weight = 84669 #+LC28
#84669
#density = 1.225
#landing = False #True: LC22, LC28
#takeoff = False #True: LC4

import sys
import os

# --- Load Case Data ---
# Structure: {'LC_Name': [load_factor, velocity, weight, density, landing, takeoff]}
LOAD_CASES = {
    "LC4_POS": [ 3.0, 60.31, 140000, 1.225, False, True],
    "LC4_NEG": [-1.5, 60.31, 140000, 1.225, False, True],
    
    "LC9_POS": [ 4.0, 107.91, 140000, 1.225, False, False],
    "LC9_NEG": [-1.5, 107.91, 140000, 1.225, False, False],  # <-- Note: Changed to -4 for your test
    
    "LC15_POS": [ 4.0, 86.52, 84669, 1.225, False, False],
    "LC15_NEG": [-1.5, 86.52, 84669, 1.225, False, False],
    
    "LC22_POS": [ 3.0, 86.3, 140000, 1.225, True, False],
    "LC22_NEG": [-1.5, 86.3, 140000, 1.225, True, False],
    
    "LC28_POS": [ 3.0, 69.19, 84669, 1.225, True, False],
    "LC28_NEG": [-1.5, 69.19, 84669, 1.225, True, False],
}

# --- User Selection Function ---

def select_load_case():
    """Prompts the user to select a load case and returns the parameters."""
    
    print("\n--- Available Load Cases ---")
    
    # Print the available options in a formatted list
    keys = list(LOAD_CASES.keys())
    for i, key in enumerate(keys):
        # Format the description for display
        n_val = LOAD_CASES[key][0]
        v_val = LOAD_CASES[key][1]
        w_val = LOAD_CASES[key][2]
        print(f"[{i+1:>2}]: {key:<10} | n={n_val:<5} | V={v_val:<7.2f} m/s | W={w_val} N")
    
    while True:
        try:
            # Prompt the user for a selection number
            selection = input("\nEnter the number of the Load Case to run (or 'q' to quit): ")
            
            if selection.lower() == 'q':
                sys.exit(0)
            
            # Convert selection to a list index
            index = int(selection) - 1
            if 0 <= index < len(keys):
                selected_key = keys[index]
                break
            else:
                print("Invalid selection. Please enter a number from the list.")
        except ValueError:
            print("Invalid input. Please enter a number or 'q'.")

    # Retrieve the selected parameters
    params = LOAD_CASES[selected_key]
    
    print(f"\nRunning Load Case: {selected_key}...")
    
    # Map the list of parameters to descriptive variable names
    (load_factor, velocity, weight, density, landing, takeoff) = params
    
    # Print the chosen values for verification
    print(f"  Load Factor (n): {load_factor}")
    print(f"  Velocity (V):    {velocity} m/s")
    print(f"  Weight (W):      {weight} N")
    print(f"  Landing Flaps:   {landing}")
    print(f"  Takeoff Flaps:   {takeoff}")
    
    return load_factor, velocity, weight, density, landing, takeoff


# --- Execute Selection and Assign Global Variables ---

# This function is called only once when the module is first imported.
try:
    (load_factor, velocity, weight, density, landing, takeoff) = select_load_case()
except Exception as e:
    # Handle the case where the script exits (e.g., user types 'q')
    print(f"Load case selection failed or exited: {e}")
    # Set safe defaults or simply let the program exit gracefully
    sys.exit(1)

# Now, any other script importing 'conditions' will have these variables defined:
# c.load_factor, c.velocity, c.weight, c.density, c.landing, c.takeoff