"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import csv
import gzip
import shutil
import requests

from _defs import *


def setup():
    create_data_dirs()
    load_pdb_report_data()
    cleanup_all_pdb_redo_dirs()


def create_data_dirs():
    print('Checking data drectories...')
    for path in (DATA_DIR, PDB_REDO_DATA_DIR, PDB_OUTPUT_DIR, MOLPROBITY_OUTPUT_DIR, PERCENTILES_OUTPUT_DIR, ROTAMER_OUTPUT_DIR, TESTING_OUTPUT_DIR, TIMING_OUTPUT_DIR):
        if not os.path.isdir(path):
            os.mkdir(path)
    print('Done.')


def get_available_pdb_ids():
    print('Getting available PDB IDs...')
    pdb_ids = [ ]
    if os.path.isdir(PDB_REDO_DATA_DIR):
        for pdb_id in os.listdir(PDB_REDO_DATA_DIR):
            if os.path.isdir(os.path.join(PDB_REDO_DATA_DIR, pdb_id)):
                gz_files = [ fn for fn in os.listdir(os.path.join(PDB_REDO_DATA_DIR, pdb_id)) if fn[-3:] == '.gz' ]
                if len(gz_files) == 4:
                    pdb_ids.append(pdb_id)
    print('Done.')
    return pdb_ids


def get_from_pdb_redo(pdb_id, output_dir):
    pdb_id = pdb_id.lower()
    check_response = requests.get(PDB_REDO_URL + pdb_id)
    if check_response.status_code != 200:
        return False
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    for suffix in PDB_REDO_SUFFIXES:
        response = requests.get(PDB_REDO_URL + pdb_id + '/' + pdb_id + suffix)
        with open(os.path.join(output_dir, pdb_id + suffix), 'wb') as outfile:
            outfile.write(response.content)
    return True


def load_pdb_report_data(force_redownload=False):
    pdb_report_values = { }
    if force_redownload or not os.path.exists(PDB_REPORT_PATH):
        print('Downloading PDB report data...')
        response = requests.get(PDB_REPORT_URL)
        if response.status_code != 200:
            print('Failed to get PDB custom report, HTTP status', response.status_code)
            return
        with open(PDB_REPORT_PATH, 'w') as outfile:
            outfile.write(response.text)
        print('Done.')
    print('Loading PDB report data...')
    with open(PDB_REPORT_PATH, 'r') as infile:
        reader = csv.reader(infile, delimiter=',', quotechar='"')
        heading_names = [ ]
        for row_id, row in enumerate(reader):
            if row_id == 0:
                heading_names = row[1:]
                continue
            pdb_id = row[0].lower()
            data = row[1:]
            pdb_report_values[pdb_id] = { }
            for heading, value in zip(heading_names, data):
                pdb_report_values[pdb_id][heading] = value
    print('Done.')
    return pdb_report_values


def decompress_pdb_redo_dir(pdb_id, suffixes=None):
    for i, suffix in enumerate(PDB_REDO_SUFFIXES):
        if suffixes is not None and i not in suffixes:
            continue
        with gzip.open(os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + suffix + '.gz'), 'rb') as infile:
            with open(os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + suffix), 'wb') as outfile:
                shutil.copyfileobj(infile, outfile)


def cleanup_pdb_redo_dir(pdb_id):
    for suffix in PDB_REDO_SUFFIXES:
        if os.path.exists(os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + suffix)):
            os.remove(os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + suffix))


def cleanup_all_pdb_redo_dirs():
    print('Cleaning up PDB REDO data directory...')
    if not os.path.isdir(PDB_REDO_DATA_DIR):
        return
    for pdb_id in os.listdir(PDB_REDO_DATA_DIR):
        if not os.path.isdir(os.path.join(PDB_REDO_DATA_DIR, pdb_id)):
            continue
        for filename in os.listdir(os.path.join(PDB_REDO_DATA_DIR, pdb_id)):
            if filename[-3:] != '.gz':
                print('Removing: ' + str(filename))
                os.remove(os.path.join(PDB_REDO_DATA_DIR, pdb_id, filename))
    print('Done.')
