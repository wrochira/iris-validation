"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os

from _defs import *
from common import load_pdb_report_data


if __name__ =='__main__':
    print('Creating data drectories...')
    for path in (DATA_DIR, PDB_REDO_DATA_DIR, PDB_OUTPUT_DIR, MOLPROBITY_OUTPUT_DIR, PERCENTILES_OUTPUT_DIR, ROTAMER_OUTPUT_DIR, TESTING_OUTPUT_DIR, TIMING_OUTPUT_DIR):
        if not os.path.isdir(path):
            os.mkdir(path)
    print('Done.')
    load_pdb_report_data()
