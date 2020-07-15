"""
Copyright 2020 William Rochira at York Structural Biology Laboratory

Included rotamer data originally from the Richardson lab
(https://github.com/rlabduke/reference_data)
"""

import os
import gzip
import pickle
import itertools

from iris_validation.utils import product
from iris_validation import SC_INCOMPLETE_STRING


LIBRARY_PATH = os.path.join(os.path.dirname(__file__), 'data', 'library.gz')
CENTRAL_VALUES_PATH = os.path.join(os.path.dirname(__file__), 'data', 'central_values.csv')

ROTAMER_LIBRARY_DATA = None
ROTAMER_CENTRAL_VALUES = None


# Unpack each byte in a byte array as four 2-bit integers
# Unsurprisingly, NumPy is much faster than pure Python for this
def unpack_bytes(in_bytes):
    try:
        # NumPy mode
        import numpy as np
        masks = np.array([ 0b11000000, 0b00110000, 0b00001100, 0b00000011 ])
        shifts = np.array([ 6, 4, 2, 0 ])
        masked = np.array(in_bytes).reshape(-1, 1) & np.array(masks)
        shifted = masked >> np.array(shifts)
        unpacked = shifted.flatten().astype('int8')
    except ImportError:
        # Python-only mode
        bits = list(itertools.product(range(4), repeat=4))
        replaced = [ bits[b] for b in in_bytes ]
        unpacked = [ bits for byte in replaced for bits in byte ]
    return unpacked


def _load_library():
    global ROTAMER_LIBRARY_DATA
    with gzip.open(LIBRARY_PATH, 'rb') as infile:
        dim_offsets, dim_bin_ranges, dim_bin_widths, dim_num_options, compressed_byte_arrays = pickle.load(infile)
    classifications = { }
    for code, compressed in compressed_byte_arrays.items():
        compressed = bytearray(compressed) # Need to cast from str in Python 2
        classifications[code] = unpack_bytes(compressed)
    ROTAMER_LIBRARY_DATA = (dim_offsets, dim_bin_ranges, dim_bin_widths, dim_num_options, classifications)


def _load_central_values():
    global ROTAMER_CENTRAL_VALUES
    central_values = { }
    with open(CENTRAL_VALUES_PATH, 'r') as infile:
        infile.readline() # Skip header line
        for line in infile.readlines():
            splitline = line.strip().split(',')
            code = splitline[0]
            rot_name = splitline[1]
            chi_means = [ float(x) for x in splitline[2:6] if x != 'None' ]
            chi_sdevs = [ float(x) for x in splitline[6:10] if x != 'None' ]
            if code not in central_values:
                central_values[code] = [ ]
            central_values[code].append((rot_name, chi_means, chi_sdevs))
    ROTAMER_CENTRAL_VALUES = central_values


def _cv_sqdiff_scores(code, chis):
    rotamer_scores = { }
    chis = [ x for x in chis if x is not None ]
    for rot_name, chi_means, chi_sdevs in ROTAMER_CENTRAL_VALUES[code]:
        chis = chis[:len(chi_means)]
        sqdiffs = [ ]
        for i in range(len(chis)):
            deltas = (chis[i]-chi_means[i], chis[i]-chi_means[i]+360)
            best_delta = deltas[int(deltas[0]>deltas[1])]
            z_score = best_delta / chi_sdevs[i]
            sqdiff = z_score**2
            sqdiffs.append(sqdiff)
        score = (sum(sqdiffs) / len(sqdiffs))**0.5
        rotamer_scores[rot_name] = score
    return rotamer_scores


def get_classification(code, chis):
    if  SC_INCOMPLETE_STRING in chis:
        return None
    if ROTAMER_LIBRARY_DATA is None:
        _load_library()
    if ROTAMER_LIBRARY_DATA is None:
        print('ERROR: failed to load rotamer library data')
        return None
    dim_offsets, dim_bin_ranges, dim_bin_widths, dim_num_options, classifications = ROTAMER_LIBRARY_DATA
    if code not in dim_offsets.keys():
        return None
    closest_values = [ ]
    chis = tuple([ x for x in chis if x is not None ][:len(dim_offsets[code])])
    for dimension, chi in enumerate(chis):
        dim_width = dim_bin_ranges[code][dimension][1] - dim_bin_ranges[code][dimension][0]
        if chi <= dim_bin_ranges[code][dimension][0]:
            chi += dim_width
        if chi >= dim_bin_ranges[code][dimension][1]:
            chi -= dim_width
        multiple = round((chi - dim_offsets[code][dimension]) / float(dim_bin_widths[code][dimension]))
        closest_value = dim_offsets[code][dimension] + multiple * dim_bin_widths[code][dimension]
        closest_values.append(closest_value)
    closest_values = tuple(closest_values)
    index = 0
    for dimension, chi in enumerate(closest_values):
        dim_offest = dim_offsets[code][dimension]
        dim_bin_width = dim_bin_widths[code][dimension]
        index += int((chi - dim_offest) / dim_bin_width * product(dim_num_options[code][dimension+1:]))
    return classifications[code][index]


def get_cv_score(code, chis):
    if ROTAMER_CENTRAL_VALUES is None:
        _load_central_values()
    if code not in ROTAMER_CENTRAL_VALUES.keys() or SC_INCOMPLETE_STRING in chis:
        return None
    scores = _cv_sqdiff_scores(code, chis)
    best_score = min(scores.values())
    return best_score
