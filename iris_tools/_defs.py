"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os

PDB_REPORT_URL = 'https://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*&customReportColumns=structureId,atomSiteCount,structureAuthor,classification,depositionDate,experimentalTechnique,macromoleculeType,ndbId,pdbDoi,releaseDate,residueCount,resolution,revisionDate,structureMolecularWeight,structureTitle&service=wsfile&format=csv'
PDB_REPORT_PATH = './data/pdb/pdb_report.csv'

PDB_REDO_URL = 'https://pdb-redo.eu/db/'
PDB_REDO_RECORD_PATH = './data/pdb/pdb_redo_record.csv'
PDB_REDO_DATA_DIR = os.path.join(os.path.expanduser('~'), 'pdb_redo_files')
PDB_REDO_SUFFIXES = ('_0cyc.pdb', '_0cyc.mtz', '_final.pdb', '_final.mtz')

DATA_DIR = './data/'
PDB_OUTPUT_DIR = './data/pdb/'
MOLPROBITY_OUTPUT_DIR = './data/molprobity/'
PERCENTILES_OUTPUT_DIR = './data/percentiles/'
ROTAMER_OUTPUT_DIR = './data/rotamer/'
TESTING_OUTPUT_DIR = './data/test_results/'
TIMING_OUTPUT_DIR = './data/timing_results/'

RESOLUTION_BIN_PERCENTILES = (10, 20, 30, 40, 50, 60, 70, 80, 90)
