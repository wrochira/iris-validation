"""
Copyright 2020 William Rochira at York Structural Biology Laboratory

- Generates an example Iris reoprt for a given PDB code using model and reflections
  data from PDB-REDO
"""

import os
import sys
import time

from iris_validation import generate_report
from iris_tools.common import get_from_pdb_redo


INPUT_DIR = 'example_input'
OUTPUT_DIR_PREFIX = 'example_report'

if __name__ == '__main__':
    pdb_id = raw_input('Enter PDB code: ')

    latest_model_path = os.path.join(INPUT_DIR, pdb_id + '_final.pdb')
    previous_model_path = os.path.join(INPUT_DIR, pdb_id + '_0cyc.pdb')
    latest_reflections_path = os.path.join(INPUT_DIR, pdb_id + '_final.mtz')
    previous_reflections_path = os.path.join(INPUT_DIR, pdb_id + '_0cyc.mtz')

    all_files_present = True
    for path in (latest_model_path, previous_model_path, latest_reflections_path, previous_reflections_path):
        if not os.path.exists(path):
            all_files_present = False
            continue

    if not all_files_present:
        if not get_from_pdb_redo(pdb_id, INPUT_DIR):
            print('Failed to get code "' + pdb_id + '" from PDB-REDO.')
            exit(1)

    generate_report(latest_model_path,
                    previous_model_path,
                    latest_reflections_path,
                    previous_reflections_path,
                    OUTPUT_DIR_PREFIX + '_' + pdb_id + '/',
                    mode='')
