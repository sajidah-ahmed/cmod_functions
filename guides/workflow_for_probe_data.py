"""

This is code in an example on how to stay organised!
Do not run this code until you have edited the code to save the files in the desired directory.

In this code I am saving mirror-Langmuir probe data.

The following shots are dwelling. You can also wish to save the rho data, which is the probe position in flux coordinates.
For the sake of this example, we will just save the time series.


"""


# Relevant imports

import numpy as np
import cmod_functions as cmod_functions


# ATTENTION: Change this path!
# Set up the folder where you wish to save the data
folder_name = f"/home/saahmed/raw_asp_mlp/"

# Specify the shots of interest
shot_list = [1160616007, 1160616018]

# Make list of probe pins
probe_pin_list = [0, 1, 2, 3]

# Make list of variable names - interrogate the available ones in the repository before you do this!
variable_list = ["Is", "ne", "Te"]


for shot in shot_list:

    for probe_pin in probe_pin_list:

        for variable in variable_list:

            # Get raw data

            time, raw_data = cmod_functions.get_raw_asp_mlp_data(
                shot, probe_pin, variable
            )

            # Set up a filename that makes sense. Here we are saving the individual parameters from each pin
            # If you printed this, for example, it'll look like: asp_mlp_1160616007_Is_pin_1.npz
            # This means that the file has ion saturation data from probe pin 1 for shot 1160616007

            file_name = f"asp_mlp_{shot}_{variable}_pin_{probe}.npz"

            np.savez(folder_name + filename, time=time, raw_data=raw_data)

            # It is extremely helpful to know whether the code has run and what it has run
            # Use print statements!

            print(f"shot {shot}, {variable}, pin {probe_pin} is complete!")
