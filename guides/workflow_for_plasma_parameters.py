"""

This is code in an example on how to stay organised!
Do not run this code until you have edited the code to save the files in the desired directory.

In this code I am saving plasma parameters for the shots of interest!


"""


# Relevant imports

import numpy as np
import cmod_functions as cmod_functions


# ATTENTION: Change this path!
# Set up the folder where you wish to save the data
folder_name = f"/home/saahmed/plasma_parameters/"

# Specify the shots of interest
shot_list = [1160616007, 1160616018]


for shot in shot_list:

    # Here, we are saving the traces of the plasma parameters for an entire discharge!

    (
        line_averaged_density_time,
        line_averaged_density,
    ) = cmod_functions.get_line_averaged_density(shot)

    plasma_current_time, plasma_current = cmod_functions.get_plasma_current(shot)

    (
        toroidal_magnetic_field_time,
        toroidal_magnetic_field,
    ) = cmod_functions.get_toroidal_magnetic_field(shot)

    file_name = f"plasma_parameters_{shot}.npz"

    np.savez(
        folder_name + file_name,
        line_averaged_density_time=line_averaged_density_time,
        line_averaged_density=line_averaged_density,
        plasma_current_time=plasma_current_time,
        plasma_current=plasma_current,
        toroidal_magnetic_field_time=toroidal_magnetic_field_time,
        toroidal_magnetic_field=toroidal_magnetic_field,
    )

    # Print a statement to say that the plasma parameters have been saved for this shot

    print(f"shot {shot} plasma parameters done!")
