"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import gzip
import time
import pickle
from multiprocessing import Process, Queue

import numpy as np
from iris_validation.metrics import generate_metrics_model
from iris_validation import METRIC_NAMES, RESOLUTION_BIN_NAMES

from _defs import PDB_REDO_DATA_DIR, PERCENTILES_OUTPUT_DIR, RESOLUTION_BIN_PERCENTILES
from common import setup, get_available_pdb_ids, load_pdb_report_data, decompress_pdb_redo_dir, cleanup_pdb_redo_dir


NUM_WORKERS = 16
YEAR_THRESHOLD = 2010
NUM_RESULTS_NEEDED = 100000000
OUTPUT_PERCENTILES = tuple([ i+1 for i in range(99) ])
EXPORT_SERIALISED = True

PDB_IDS = [ ]
PDB_REPORT_DATA = { }
RESOLUTION_BINS = { }
METRIC_VALUES_ALL = { }
METRIC_VALUES_BINNED = { }
METRIC_PERCENTILE_VALUES = { }
NUM_MODELS_ANALYSED = 0


def get_resolution_bin_id(resolution):
    bin_id = 9
    for i, percentile in enumerate(sorted(RESOLUTION_BINS.keys())):
        percentile_resolution = RESOLUTION_BINS[percentile]
        if resolution < percentile_resolution:
            bin_id = i
            break
    return bin_id


def worker(in_queue, out_queue):
    while True:
        if in_queue.empty():
            exit()
        pdb_id = in_queue.get()
        metric_values = { }
        for metric_name in METRIC_NAMES:
            metric_values[metric_name] = [ ]
        decompress_pdb_redo_dir(pdb_id)
        model_path = os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + '_0cyc.pdb')
        reflections_path = os.path.join(PDB_REDO_DATA_DIR, pdb_id, pdb_id + '_0cyc.mtz')
        try:
            metrics_model = generate_metrics_model(model_path, reflections_path)
        except:
            continue
        cleanup_pdb_redo_dir(pdb_id)
        for chain in metrics_model:
            for residue in chain:
                # Skip non amino acid residues
                if not residue.is_aa:
                    continue
                metric_values['Ramachandran Score'].append(residue.ramachandran_score)
                metric_values['Rotamer Score'].append(residue.rotamer_score)
                metric_values['Avg B-factor'].append(residue.avg_b_factor)
                metric_values['Max B-factor'].append(residue.max_b_factor)
                metric_values['Std B-factor'].append(residue.std_b_factor)
                metric_values['Residue Fit'].append(residue.fit_score)
                metric_values['Mainchain Fit'].append(residue.mainchain_fit_score)
                metric_values['Sidechain Fit'].append(residue.sidechain_fit_score)
        # Discard null values
        for metric_name in metric_values.keys():
            metric_values[metric_name] = [ x for x in metric_values[metric_name] if x is not None ]
        out_queue.put((pdb_id, metric_values))


def generate_percentiles_data():
    global RESOLUTION_BINS
    global METRIC_VALUES_ALL
    global METRIC_VALUES_BINNED
    global METRIC_PERCENTILE_VALUES
    global NUM_MODELS_ANALYSED
    # Add all the avilable PDB IDs to the input queue
    in_queue, out_queue = Queue(), Queue()
    pdb_resolutions = { }
    for pdb_id in PDB_IDS:
        pdb_data = PDB_REPORT_DATA[pdb_id]
        # Skip if deposited prior to the year threshold
        year = int(pdb_data['depositionDate'].split('-')[0])
        if YEAR_THRESHOLD is not None and year < YEAR_THRESHOLD:
            continue
        # Record the resolution for later classification, or skip the model if the resolution is null
        try:
            pdb_resolutions[pdb_id] = float(pdb_data['resolution'])
        except:
            continue
        in_queue.put(pdb_id)
    print('*** Number of structures: ' + str(len(pdb_resolutions.keys())))
    # Calculate the resolution bins based on all the recorded resolutions
    resolution_bin_values = np.percentile(pdb_resolutions.values(), RESOLUTION_BIN_PERCENTILES)
    for bin_percentile, bin_value in zip(RESOLUTION_BIN_PERCENTILES, resolution_bin_values):
        RESOLUTION_BINS[bin_percentile] = bin_value
    print(RESOLUTION_BINS)
    # Initialise METRIC_VALUES_BINNED as a 2D dictionary indexed by metric name and then resolution bin ID
    for metric_name in METRIC_NAMES:
        METRIC_VALUES_BINNED[metric_name] = { }
        for resolution_bin_id in range(len(RESOLUTION_BIN_PERCENTILES)+1):
            METRIC_VALUES_BINNED[metric_name][resolution_bin_id] = [ ]
    print('Generating metrics...')
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
        num_results = [ sum([ len(y) for y in x.values() ]) for x in METRIC_VALUES_BINNED.values() ]
        print('*** Number of results: ' + str(num_results))
        # Stop if desired number of results have been calculated
        if NUM_RESULTS_NEEDED is not None and min(num_results) >= NUM_RESULTS_NEEDED:
            for p in processes:
                p.terminate()
            break
        # Stop if all IDs have been processed
        if in_queue.empty() and out_queue.empty():
            break
        # Collect results from the output queue
        while not out_queue.empty():
            pdb_id, results = out_queue.get()
            resolution = pdb_resolutions[pdb_id]
            resolution_bin_id = get_resolution_bin_id(resolution)
            for metric_name, values in results.items():
                METRIC_VALUES_BINNED[metric_name][resolution_bin_id] += values
            NUM_MODELS_ANALYSED += 1
        time.sleep(1)
    print('Calculating percentiles...')
    for metric_name in METRIC_NAMES:
        METRIC_VALUES_ALL[metric_name] = [ ]
        METRIC_PERCENTILE_VALUES[metric_name] = { }
        for resolution_bin_id in range(len(RESOLUTION_BIN_PERCENTILES)+1):
            resolution_metric_values = METRIC_VALUES_BINNED[metric_name][resolution_bin_id]
            METRIC_VALUES_ALL[metric_name] += resolution_metric_values
            if len(resolution_metric_values) == 0:
                print('ERROR: Empty bin! ' + str(metric_name + ', bin ID ' + str(resolution_bin_id)))
                continue
            percentile_values = np.percentile(resolution_metric_values, OUTPUT_PERCENTILES)
            METRIC_PERCENTILE_VALUES[metric_name][resolution_bin_id] = percentile_values
        METRIC_PERCENTILE_VALUES[metric_name]['All'] = np.percentile(METRIC_VALUES_ALL[metric_name], OUTPUT_PERCENTILES)
    print('Done.')

def export_percentiles_data():
    if not os.path.isdir(PERCENTILES_OUTPUT_DIR):
        os.mkdir(PERCENTILES_OUTPUT_DIR)
    print('Exporting resolution bin thresholds...')
    with open(os.path.join(PERCENTILES_OUTPUT_DIR, 'resolution_bins.csv'), 'w') as outfile:
        outfile.write('Bin Name,Resolution\n')
        for percentile in RESOLUTION_BIN_PERCENTILES:
            threshold = RESOLUTION_BINS[percentile]
            outfile.write(str(percentile) + ',' + str(round(threshold, 3)) + '\n')
    print('Exporting percentiles data...')
    with open(os.path.join(PERCENTILES_OUTPUT_DIR, 'percentiles_data.csv'), 'w') as outfile:
        outfile.write(','.join([ 'Resolution Bin', 'Percentile' ] + list(METRIC_NAMES)) + '\n')
        for resolution_bin_id in range(len(RESOLUTION_BIN_PERCENTILES)+1):
            for percentile in OUTPUT_PERCENTILES:
                line_values = [ RESOLUTION_BIN_NAMES[resolution_bin_id], percentile ]
                for metric_name in METRIC_NAMES:
                    value = METRIC_PERCENTILE_VALUES[metric_name][resolution_bin_id][percentile-1]
                    line_values.append(round(value, 5))
                outfile.write(','.join([ str(x) for x in line_values ]) + '\n')
        for percentile in OUTPUT_PERCENTILES:
            line_values = [ 'All', percentile ]
            for metric_name in METRIC_NAMES:
                value = METRIC_PERCENTILE_VALUES[metric_name]['All'][percentile-1]
                line_values.append(round(value, 5))
            outfile.write(','.join([ str(x) for x in line_values ]) + '\n')
    print('Exporting sample sizes...')
    with open(os.path.join(PERCENTILES_OUTPUT_DIR, 'sample_sizes.csv'), 'w') as outfile:
        outfile.write('Metric name,Sample size (residues),Sample size (models)\n')
        for metric_name in METRIC_NAMES:
            num_results = len(METRIC_VALUES_ALL[metric_name])
            outfile.write(metric_name + ',' + str(num_results) + ',' + str(NUM_MODELS_ANALYSED) + '\n')
    if EXPORT_SERIALISED:
        print('Exporting serialised data...')
        with gzip.open(os.path.join(PERCENTILES_OUTPUT_DIR, 'serialised_metrics.gz'), 'wb') as outfile:
            pickle.dump(METRIC_VALUES_BINNED, outfile)
    print('Done.')


if __name__ == '__main__':
    setup()
    PDB_IDS = get_available_pdb_ids()
    PDB_REPORT_DATA = load_pdb_report_data()
    generate_percentiles_data()
    export_percentiles_data()
