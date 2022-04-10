# cmod_functions
Functions for extracting data from Alcator C-Mod while on the MFE cluster.

# Prerequisites
eqtools: https://github.com/PSFCPlasmaTools/eqtools \
Please be aware that eqtools may not be fully Python3 compatible.
````
cd eqtools/eqtools
vim core.py
# Replace line number 45: import trispline
# With: from . import trispline
# Save the file.
````

MDSplus: https://www.mdsplus.org/

Xarray: https://docs.xarray.dev/en/stable/getting-started-guide/installing.html


# Installation
```
git clone https://github.com/sajidah-ahmed/cmod_functions.git
cd cmod_functions
pip install -e .
```
