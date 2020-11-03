# Iris-Validation
A Python package for interactive all-in-one graphical validation of 3D protein model iterations.

## metrics
Calculates comprehensive per-residue metrics from a model (and reflections) file.

## interface
Generates the HTML validation report.

## utils
A number of functions that are utilised extensively within the library and are also useful on their own.


# Iris-Tools
Companion scripts for the main module that include tests, timings, and the scripts used to generate data for the metrics module.


# Installation
The package can now be installed via pip, with the command `pip install iris-validation`.
For now, the best way to install and run Iris is to use the Python environment that comes with CCP4, by running `ccp4-python -m pip install iris-validation`
You may first need to install pip, by downloading the install script from https://bootstrap.pypa.io/get-pip.py and running `ccp4-python get-pip.py`
This is because Clipper-Python, upon which Iris depends, is currently available in a number of different forms, most of which are broken in some way or another. As a consequece, only the Python 2 version of Iris is available. This will be remedied soon.
