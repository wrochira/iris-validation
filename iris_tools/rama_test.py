import os
import gzip
import time
import pickle
import shutil
import itertools
from multiprocessing import Process, Queue

from iotbx import file_reader
from cStringIO import StringIO
from iris_validation import utils
from mmtbx.validation import ramalyze

from _defs import PDB_REDO_DATA_DIR, TESTING_OUTPUT_DIR
from common import setup, get_available_pdb_ids, load_pdb_report_data, decompress_pdb_redo_dir, cleanup_pdb_redo_dir, cleanup_all_pdb_redo_dirs


NUM_WORKERS = 6
YEAR_THRESHOLD = None
CLASSIFICATIONS = { 'O' : 1, 'A' : 2, 'F' : 3 }
CLASSIFICATION_PAIRS = list(itertools.product((1, 2, 3), (1, 2, 3)))

PDB_IDS = [ ]
PDB_REPORT_DATA = { }
CLASS_RESULTS = { }
SCORE_RESULTS = { }


def worker(in_queue, out_queue):
    while True:
        if in_queue.empty():
            exit()
        pdb_id = in_queue.get()
        decompress_pdb_redo_dir(pdb_id, suffixes={0})
        pdb_path = os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + '_0cyc.pdb')
        pdb_in = file_reader.any_file(file_name=pdb_path)
        hierarchy = pdb_in.file_object.hierarchy
        rama_analysis = ramalyze.ramalyze(pdb_hierarchy=hierarchy, outliers_only=False)
        out = StringIO()
        rama_analysis.show_old_output(out=out, verbose=False)
        output = out.getvalue()

        class_results = { }
        score_results = { }
        for pair in CLASSIFICATION_PAIRS:
            class_results[pair] = 0
        for clf in (1, 2, 3):
            score_results[clf] = [ ]

        for line in output.split('\n'):
            if len(line) == 0:
                continue
            chain_id = line[:2].strip()
            seqnum = int(line[2:6].strip())
            splitline = [ x.strip() for x in line[6:].split(':') ]
            code = splitline[0].strip()
            if code == 'MSE':
                code = 'MET'
            phi, psi = [ float(x) for x in splitline[2:4] if len(x) > 0 ]
            iris_score = utils.calculate_ramachandran_score(None, code, phi, psi)
            iris_class_int = 1 if iris_score < 0.002 else 2 if iris_score < 0.02 else 3
            if code not in utils.ONE_LETTER_CODES.values():
                continue
            if iris_class_int == None:
                continue
            mp_class_char = splitline[-2][0].upper()
            mp_class_int = CLASSIFICATIONS[mp_class_char]
            class_results[(iris_class_int, mp_class_int)] += 1
            score_results[mp_class_int].append(iris_score)
        cleanup_pdb_redo_dir(pdb_id)
        out_queue.put((pdb_id, class_results, score_results))


def ramalyze_all():
    global CLASS_RESULTS
    for pair in CLASSIFICATION_PAIRS:
        CLASS_RESULTS[pair] = 0
    for clf in (1, 2, 3):
        SCORE_RESULTS[clf] = [ ]
    # Add all the avilable PDB IDs to the input queue
    in_queue, out_queue = Queue(), Queue()
    for pdb_id in PDB_IDS:
        pdb_data = PDB_REPORT_DATA[pdb_id]
        # Skip if deposited prior to the year threshold
        year = int(pdb_data['depositionDate'].split('-')[0])
        if YEAR_THRESHOLD is not None and year < YEAR_THRESHOLD:
            continue
        in_queue.put(pdb_id)
    print('Ramalyzing...')
    processes = [ ]
    num_models_analysed = 0
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
            pdb_id, class_results, score_results = out_queue.get()
            for pair in CLASSIFICATION_PAIRS:
                CLASS_RESULTS[pair] += class_results[pair]
            for clf in (1, 2, 3):
                SCORE_RESULTS[clf] += score_results[clf]
            num_models_analysed += 1
        print('*** Models analysed: ' + str(num_models_analysed))
        time.sleep(1)
    print('Exporting confusion matrix...')
    with open(os.path.join(TESTING_OUTPUT_DIR, 'rama_clf.csv'), 'w') as outfile:
        outfile.write(',,Iris,,\n')
        outfile.write(',,Outlier,Allowed,Favoured\n')
        outfile.write('Molprobity,Outlier,' + ','.join([ str(CLASS_RESULTS[(i+1, 1)]) for i in range(3) ]) + '\n')
        outfile.write(',Allowed,' + ','.join([ str(CLASS_RESULTS[(i+1, 2)]) for i in range(3) ]) + '\n')
        outfile.write(',Favoured,' + ','.join([ str(CLASS_RESULTS[(i+1, 3)]) for i in range(3) ]) + '\n')
        outfile.write('\n')
        outfile.write('Num models,' + str(num_models_analysed))
    print('Exporting score analyses...')
    import numpy as np
    with open(os.path.join(TESTING_OUTPUT_DIR, 'rama_scores_summary.csv'), 'w') as outfile:
        for clf in (1, 2, 3):
            mean = np.mean(SCORE_RESULTS[clf])
            std = np.std(SCORE_RESULTS[clf])
            outfile.write(str(clf) + ',' + str(mean) + ',' + str(std) + '\n')
    with open(os.path.join(TESTING_OUTPUT_DIR, 'rama_scores_all.csv'), 'w') as outfile:
        outfile.write('1,2,3\n')
        results_counts = [ len(x) for x in SCORE_RESULTS.values() ]
        for i in range(max(results_counts)):
            outline = ''
            for clf in (1, 2, 3):
                if len(SCORE_RESULTS[clf]) >= i+1:
                    outline += str(SCORE_RESULTS[clf][i])
                outline += ','
            outfile.write(outline + '\n')
    print('Done.')


if __name__ == '__main__':
    setup()
    PDB_IDS = get_available_pdb_ids()
    PDB_REPORT_DATA = load_pdb_report_data()
    ramalyze_all()
    cleanup_all_pdb_redo_dirs()
