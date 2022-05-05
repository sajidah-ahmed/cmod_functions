# cmod_functions
Functions for extracting data from Alcator C-Mod while on the MFE cluster. \
Guides show how to analyse some of the data or use the functions if you're doing analysis on the cluster or on your local machine.
* How to use the xarray function to generate APD data.
* Plot profiles with scanning mirror-Langmuir probe data.
* How to download Phantom camera data and turn it into an animation

# Prerequisites
MDSplus: https://www.mdsplus.org/

eqtools: https://github.com/PSFCPlasmaTools/eqtools \
Please be aware that eqtools may not be fully Python3-compatible.
````
cd eqtools/eqtools
vim core.py
# Replace line number 45: import trispline
# With: from . import trispline
# Save the file.
````
NOTE: To use the functions involving eqtools, you must have some knowledge of which shots to run it for and for which EFIT time resolution.
According to Jim Terry, there exists a relatively poor resolution of the default EFIT ("ANALYSIS" with a 33x33 spatial grid, and 20 ms time resolution). If you want a higher resolution EFIT (129x129 grid with a time resolution as small as 0.5 ms), you must inform someone about this. 

# Installation
```
git clone https://github.com/sajidah-ahmed/cmod_functions.git
cd cmod_functions
pip install -e .
```
