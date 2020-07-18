import os
import gzip
import time
import pickle
import shutil
import itertools
from math import pi
from multiprocessing import Process, Queue

import numpy as np
from iotbx import file_reader
from cStringIO import StringIO
import matplotlib.pyplot as plt
from iris_validation import utils
from mmtbx.validation import ramalyze

from _defs import PDB_REDO_DATA_DIR, TESTING_OUTPUT_DIR
from common import setup, get_available_pdb_ids, load_pdb_report_data, decompress_pdb_redo_dir, cleanup_pdb_redo_dir, cleanup_all_pdb_redo_dirs


NUM_WORKERS = 8
YEAR_THRESHOLD = None
THRESHOLDS = (0.01, 0.0005) # Clipper default
#THRESHOLDS = (0.02, 0.002) # Concordant with Coot
CLASSIFICATIONS = ('Outlier', 'Allowed', 'Favored')
CLASSIFICATION_PAIRS = list(itertools.product(CLASSIFICATIONS, repeat=2))

PDB_IDS = [ ]
PDB_REPORT_DATA = { }
CLASS_RESULTS = { }
SCORE_RESULTS = { }
NUM_MODELS_ANALYSED = 0


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
        for pair in CLASSIFICATION_PAIRS:
            class_results[pair] = 0
        score_results = { }
        for clf in CLASSIFICATIONS:
            score_results[clf] = [ ]

        for line in output.split('\n'):
            if len(line) == 0:
                continue
            chain_id = line[:2].strip()
            seqnum = int(line[2:6].strip())
            splitline = [ x.strip() for x in line[6:].split(':') ]
            code = splitline[0].strip()
            if code not in utils.THREE_LETTER_CODES[0]:
                continue
            if code == 'MSE':
                code = 'MET'
            phi, psi = [ float(x) for x in splitline[2:4] ]
            phi *= pi/180
            psi *= pi/180
            iris_score = utils.calculate_ramachandran_score(None, code, phi, psi)
            iris_class = CLASSIFICATIONS[0] if iris_score < THRESHOLDS[1] else CLASSIFICATIONS[1] if iris_score < THRESHOLDS[0] else CLASSIFICATIONS[2]
            mp_class = splitline[-2][0].upper() + splitline[-2][1:].lower()
            class_results[(iris_class, mp_class)] += 1
            score_results[mp_class].append(iris_score)
        cleanup_pdb_redo_dir(pdb_id)
        out_queue.put((pdb_id, class_results, score_results))


def ramalyze_all():
    global CLASS_RESULTS
    global NUM_MODELS_ANALYSED
    for pair in CLASSIFICATION_PAIRS:
        CLASS_RESULTS[pair] = 0
    for clf in CLASSIFICATIONS:
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
        # Collect results from the output queue
        while not out_queue.empty():
            pdb_id, class_results, score_results = out_queue.get()
            for pair in CLASSIFICATION_PAIRS:
                CLASS_RESULTS[pair] += class_results[pair]
            for clf in CLASSIFICATIONS:
                SCORE_RESULTS[clf] += score_results[clf]
            NUM_MODELS_ANALYSED += 1
        print('*** Models analysed: ' + str(NUM_MODELS_ANALYSED))
        time.sleep(1)


def export_results():
    print('Exporting confusion matrix...')
    with open(os.path.join(TESTING_OUTPUT_DIR, 'rama_clf.csv'), 'w') as outfile:
        outfile.write(',,Iris,,\n')
        outfile.write(',,' + ','.join(CLASSIFICATIONS) + '\n')
        outfile.write('Molprobity')
        for clf_molprobity in CLASSIFICATIONS:
            outfile.write(',' + clf_molprobity + ',' + ','.join([ str(CLASS_RESULTS[(clf_iris, clf_molprobity)]) for clf_iris in CLASSIFICATIONS ]) + '\n')
        outfile.write('\nNum models,' + str(NUM_MODELS_ANALYSED))
    print('Exporting score analyses...')
    with open(os.path.join(TESTING_OUTPUT_DIR, 'rama_scores_summary.csv'), 'w') as outfile:
        for clf in CLASSIFICATIONS:
            mean = np.mean(SCORE_RESULTS[clf])
            std = np.std(SCORE_RESULTS[clf])
            outfile.write(clf + ',' + str(mean) + ',' + str(std) + '\n')
    with open(os.path.join(TESTING_OUTPUT_DIR, 'rama_scores_all.csv'), 'w') as outfile:
        outfile.write(','.join(CLASSIFICATIONS) + '\n')
        results_counts = [ len(x) for x in SCORE_RESULTS.values() ]
        for i in range(max(results_counts)):
            outline = ''
            for clf in CLASSIFICATIONS:
                if len(SCORE_RESULTS[clf]) >= i+1:
                    outline += str(SCORE_RESULTS[clf][i])
                outline += ','
            outfile.write(outline + '\n')
    print('Done.')


def draw_graphs():
    print('Making histograms...')
    for clf in CLASSIFICATIONS:
        plt.hist(SCORE_RESULTS[clf], bins=100)
        plt.title('Molprobity ' + clf)
        ax = plt.gca()
        ax.set_ylabel('Count')
        ax.set_xlabel('Clipper Score')
        ax.set_yscale('log')
        plt.savefig(os.path.join(TESTING_OUTPUT_DIR, 'rama_hist_' + clf.lower() + '.png'), dpi=600)
        plt.close()


if __name__ == '__main__':
    setup()
    PDB_IDS = get_available_pdb_ids()
    PDB_REPORT_DATA = load_pdb_report_data()
    ramalyze_all()
    export_results()
    draw_graphs()
    cleanup_all_pdb_redo_dirs()
