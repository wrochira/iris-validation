"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import gzip
import pickle

import matplotlib.pyplot as plt
from iris_validation import METRIC_NAMES

from _defs import PERCENTILES_OUTPUT_DIR, RESOLUTION_BIN_PERCENTILES


METRIC_VALUES_ALL = { }
METRIC_VALUES_BINNED = { }

PRESET_X_MAXIMA = { 'Ramachandran Score' : 2.00,
                    'Rotamer Score' : 17,
                    'Avg B-factor' : 500,
                    'Max B-factor' : 500,
                    'Std B-factor' : 30,
                    'Residue Fit' : 1000,
                    'Mainchain Fit' : 1000,
                    'Sidechain Fit' : 1000 }
USE_PRESETS = True
USE_X_LIMITS = True
NUM_BINS = 1000


def load_data():
    global METRIC_VALUES_ALL
    global METRIC_VALUES_BINNED

    with gzip.open(os.path.join(PERCENTILES_OUTPUT_DIR, 'serialised_metrics.gz'), 'rb') as infile:
        METRIC_VALUES_BINNED = pickle.load(infile)
    for metric_name in METRIC_NAMES:
        METRIC_VALUES_ALL[metric_name] = [ ]
        for resolution_bin_id in range(len(RESOLUTION_BIN_PERCENTILES)+1):
            resolution_metric_values = METRIC_VALUES_BINNED[metric_name][resolution_bin_id]
            METRIC_VALUES_ALL[metric_name] += resolution_metric_values


def draw_graphs():
    if not os.path.isdir(os.path.join(PERCENTILES_OUTPUT_DIR, 'graphs')):
        os.mkdir(os.path.join(PERCENTILES_OUTPUT_DIR, 'graphs'))

    # Get maximum values to set axis limits
    if USE_PRESETS:
        x_maxima = PRESET_X_MAXIMA
    else:
        x_maxima = { }
        for metric_name in METRIC_NAMES:
            x_maxima[metric_name] = max(METRIC_VALUES_BINNED[metric_name][resolution_bin_id][-1] for resolution_bin_id in range(len(RESOLUTION_BIN_PERCENTILES)+1))

    # Generate and export histograms
    print('Generating overview histograms...')
    for metric_name in METRIC_NAMES:
        print('*** ' + metric_name)
        plt.hist(METRIC_VALUES_ALL[metric_name], bins=NUM_BINS)
        plt.title(metric_name)
        ax = plt.gca()
        ax.set_ylabel('Count')
        ax.set_xlabel('Value')
        if USE_X_LIMITS:
            ax.set_xlim((0, x_maxima[metric_name]))
        plt.savefig(os.path.join(PERCENTILES_OUTPUT_DIR, 'graphs', metric_name + '_hist_all' + '.png'), dpi=600)
        plt.close()
    print('Generating binned histograms...')
    for metric_name in METRIC_NAMES:
        print('*** ' + metric_name)
        for resolution_bin_id in range(len(RESOLUTION_BIN_PERCENTILES)+1):
            values = METRIC_VALUES_BINNED[metric_name][resolution_bin_id]
            plt.hist(values, bins=NUM_BINS, alpha=0.25)
        plt.title(metric_name)
        ax = plt.gca()
        ax.set_ylabel('Count')
        ax.set_xlabel('Value')
        if USE_X_LIMITS:
            ax.set_xlim((0, x_maxima[metric_name]))
        plt.savefig(os.path.join(PERCENTILES_OUTPUT_DIR, 'graphs', metric_name + '_hist_combo.png'), dpi=600)
        plt.close()
        for resolution_bin_id in range(len(RESOLUTION_BIN_PERCENTILES)+1):
            values = METRIC_VALUES_BINNED[metric_name][resolution_bin_id]
            plt.hist(values, bins=NUM_BINS)
            plt.title(metric_name)
            ax = plt.gca()
            ax.set_ylabel('Count')
            ax.set_xlabel('Value')
            if USE_X_LIMITS:
                ax.set_xlim((0, x_maxima[metric_name]))
            plt.savefig(os.path.join(PERCENTILES_OUTPUT_DIR, 'graphs', metric_name + '_hist_bin_' + str(resolution_bin_id) + '.png'), dpi=600)
            plt.close()
    print('Done.')


if __name__ == '__main__':
    load_data()
    draw_graphs()
