"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import gzip
import time
import pickle
import shutil
import itertools
from multiprocessing import Process, Queue

from iotbx import file_reader
from cStringIO import StringIO
from mmtbx.validation import rotalyze
from iris_validation.metrics import rotamer

from _defs import PDB_REDO_DATA_DIR, TESTING_OUTPUT_DIR
from common import setup, get_available_pdb_ids, load_pdb_report_data, decompress_pdb_redo_dir, cleanup_pdb_redo_dir, cleanup_all_pdb_redo_dirs


NUM_WORKERS = 16
YEAR_THRESHOLD = None
CLASSIFICATIONS = { 'O' : 1, 'A' : 2, 'F' : 3 }
CLASSIFICATION_PAIRS = list(itertools.product((0, 1, 2, 3), (1, 2, 3)))

PDB_IDS = [ ]
PDB_REPORT_DATA = { }
RESULTS = { }
NUM_MODELS_ANALYSED = 0


def worker(in_queue, out_queue):
    while True:
        if in_queue.empty():
            exit(0)
        pdb_id = in_queue.get()
        decompress_pdb_redo_dir(pdb_id, suffixes={0})
        pdb_path = os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + '_0cyc.pdb')
        pdb_in = file_reader.any_file(file_name=pdb_path)
        hierarchy = pdb_in.file_object.hierarchy
        rota_analysis = rotalyze.rotalyze(pdb_hierarchy=hierarchy, outliers_only=False)
        out = StringIO()
        rota_analysis.show_old_output(out=out, verbose=False)
        output = out.getvalue()
        results = { }
        for pair in CLASSIFICATION_PAIRS:
            results[pair] = 0
        for line in output.split('\n'):
            if len(line) == 0:
                continue
            chain_id = line[:2].strip()
            seqnum = int(line[2:6].strip())
            splitline = [ x.strip() for x in line[6:].split(':') ]
            code = splitline[0]
            chis = [ float(x) for x in splitline[3:7] if len(x) > 0 ]
            iris_class_int = rotamer.get_classification(code, chis)
            if iris_class_int == None:
                continue
            mp_class_char = splitline[-2][0].upper()
            mp_class_int = CLASSIFICATIONS[mp_class_char]
            results[(iris_class_int, mp_class_int)] += 1
        cleanup_pdb_redo_dir(pdb_id)
        out_queue.put((pdb_id, results))


def rotalyze_all():
    global RESULTS
    global NUM_MODELS_ANALYSED
    for pair in CLASSIFICATION_PAIRS:
        RESULTS[pair] = 0
    # Add all the avilable PDB IDs to the input queue
    in_queue, out_queue = Queue(), Queue()
    for pdb_id in PDB_IDS:
        pdb_data = PDB_REPORT_DATA[pdb_id]
        # Skip if deposited prior to the year threshold
        year = int(pdb_data['depositionDate'].split('-')[0])
        if YEAR_THRESHOLD is not None and year < YEAR_THRESHOLD:
            continue
        in_queue.put(pdb_id)
    print('Rotalyzing...')
    processes = [ ]
    while True:
        # Delete references to processes that have died (due to uncatchable Clipper errors)
        processes = [ p for p in processes if p.is_alive() ]
        # Spawn processes until NUM_WORKERS is reached
        if not in_queue.empty():
            while len(processes) < NUM_WORKERS:
                print('*** Spawning new process...')
                new_process = Process(target=worker, args=(in_queue, out_queue))
                new_process.start()
                processes.append(new_process)
        # Stop if all IDs have been processed
        if in_queue.empty() and out_queue.empty():
            break
        # Collect results from the output Queue
        while not out_queue.empty():
            pdb_id, results = out_queue.get()
            for pair in CLASSIFICATION_PAIRS:
                RESULTS[pair] += results[pair]
            NUM_MODELS_ANALYSED += 1
        print('*** Models analysed: ' + str(NUM_MODELS_ANALYSED))
        time.sleep(1)


def export_results():
    print('Exporting confusion matrix...')
    with open(os.path.join(TESTING_OUTPUT_DIR, 'rotamer.csv'), 'w') as outfile:
        outfile.write(',,Iris,,,\n')
        outfile.write(',,Unknown,Outlier,Allowed,Favoured\n')
        outfile.write('Molprobity,Outlier,' + ','.join([ str(RESULTS[(i, 1)]) for i in range(4) ]) + '\n')
        outfile.write(',Allowed,' + ','.join([ str(RESULTS[(i, 2)]) for i in range(4) ]) + '\n')
        outfile.write(',Favoured,' + ','.join([ str(RESULTS[(i, 3)]) for i in range(4) ]) + '\n')
        outfile.write('\n')
        outfile.write('Num models,' + str(NUM_MODELS_ANALYSED))
    print('Done.')


if __name__ == '__main__':
    setup()
    PDB_IDS = get_available_pdb_ids()
    PDB_REPORT_DATA = load_pdb_report_data()
    rotalyze_all()
    export_results()
    cleanup_all_pdb_redo_dirs()
