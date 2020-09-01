"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import csv
import gzip
import requests
from multiprocessing import Process, Queue

from common import setup, get_available_pdb_ids, load_pdb_report_data
from _defs import PDB_REDO_DATA_DIR, PDB_REDO_SUFFIXES, PDB_REDO_UR, PDB_REDO_RECORD_PATH


NUM_WORKERS = 16

PDB_REPORT_VALUES = { }
PDB_REDO_RECORD = { }


def worker(queue):
    while True:
        pdb_id = queue.get()
        structure_dir = os.path.join(PDB_REDO_DATA_DIR, pdb_id)
        if not os.path.isdir(structure_dir):
            os.mkdir(structure_dir)
        for suffix in PDB_REDO_SUFFIXES:
            response = requests.get(PDB_REDO_URL + pdb_id + '/' + pdb_id + suffix)
            with gzip.open(os.path.join(structure_dir, pdb_id + suffix + '.gz'), 'wb') as outfile:
                outfile.write(response.content)


def scrape_pdb_redo(retry_failed_ids=True):
    global PDB_REDO_RECORD
    # Make output directory if it does not exist
    if not os.path.exists(PDB_REDO_DATA_DIR):
        os.mkdir(PDB_REDO_DATA_DIR)
    # Get PDB IDs from PDB report, and alphabetise
    pdb_ids = sorted(PDB_REPORT_VALUES.keys())
    # Initialise record and load existing record from file if it exists
    for pdb_id in pdb_ids:
        PDB_REDO_RECORD[pdb_id] = 0
    if os.path.exists(PDB_REDO_RECORD_PATH):
        with open(PDB_REDO_RECORD_PATH, 'r') as infile:
            for row in csv.reader(infile, delimiter=','):
                pdb_id, status = row
                PDB_REDO_RECORD[pdb_id] = int(status)
    # Start workers
    queue = Queue()
    processes = [ Process(target=worker, args=(queue, )) for _ in range(NUM_WORKERS) ]
    for p in processes:
        p.start()
    # For each PDB ID, check if it needs to be checked, and if so, add it to the appropriate shortlist
    shortlists = [ [ ], [ ], [ ] ]
    for pdb_id in reversed(pdb_ids):
        status = PDB_REDO_RECORD[pdb_id]
        if status == 0:
            shortlists[0].append(pdb_id)
        elif status == 200:
            if [ os.path.exists(os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + suffix + '.gz')) for suffix in PDB_REDO_SUFFIXES ].count(True) < 4:
                shortlists[1].append(pdb_id)
        elif retry_failed_ids:
            shortlists[2].append(pdb_id)
    print('Unchecked IDs:', len(shortlists[0]))
    print('Successes to retry:', len(shortlists[1]))
    print('Failures to retry:', len(shortlists[2]))
    compiled_shortlist = [ pdb_id for shortlist in shortlists for pdb_id in shortlist ]
    shortlist_length = len(compiled_shortlist)
    # For each PDB ID in the shortlist, check if it exists in PDB-REDO, and if so, add it to the download queue
    for i, pdb_id in enumerate(compiled_shortlist):
        response = requests.get(PDB_REDO_URL + pdb_id)
        PDB_REDO_RECORD[pdb_id] = response.status_code
        print(pdb_id + ' '  + str(response.status_code) + ' ::: ' + str(i) + ' / ' + str(shortlist_length))
        if response.status_code == 200:
            queue.put(pdb_id)
        # Every 10 codes, write out the record file
        if i % 10 == 0 or i == shortlist_length-1:
            with open(PDB_REDO_RECORD_PATH, 'w') as outfile:
                for pdb_id in pdb_ids:
                    outfile.write(pdb_id + ',' + str(PDB_REDO_RECORD[pdb_id]) + '\n')


if __name__ == '__main__':
    setup()
    PDB_REPORT_VALUES = load_pdb_report_data()
    scrape_pdb_redo()
    print('Done.')
