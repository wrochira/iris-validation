"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import sys
import time

from iris_validation import generate_report


if __name__ == '__main__':
    sys.dont_write_bytecode = True

    latest_model_path = os.path.join('example_input', '2a0x_final.pdb')
    previous_model_path = os.path.join('example_input', '2a0x_0cyc.pdb')
    latest_reflections_path = os.path.join('example_input', '2a0x_final.mtz')
    previous_reflections_path = os.path.join('example_input', '2a0x_0cyc.mtz')

generate_report(latest_model_path,
                previous_model_path,
                latest_reflections_path,
                previous_reflections_path,
                './example_report/',
                mode='')
