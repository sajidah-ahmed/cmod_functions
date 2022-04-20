import os
from setuptools import setup

name = "cmod_functions"

with open("README.md") as f:
    long_description = f.read()

here = os.path.abspath(os.path.dirname(__file__))

setup(
    name=name,
    description="Python scripts for Alcator C-Mod data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sajidah-ahmed/cmod_functions",
    license="MiT",
    version="1.0",
    packages=["cmod_functions"],
    python_requires=">=3.0",
    install_requires=[
        "numpy>=1.18.0",
        "scipy>=1.4.0",
        "matplotlib>=3.2.0",
        "xarray>=0.16.2",
        "dask",
        "netCDF4",
        "bottleneck",
    ],
    classifiers=[
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    zip_safe=False,
)
