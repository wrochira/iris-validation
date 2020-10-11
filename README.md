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
At the moment, the Iris package is not compatible with Python 3 under Windows. This is because Clipper-Python, upon which Iris depends, is not currently available for that environment. This will hopefully be remedied soon.
