"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import gzip
import pickle
import zipfile
import requests

import numpy as np

from common import setup
from _defs import ROTAMER_OUTPUT_DIR


REFERENCE_DATA_URL = 'https://github.com/rlabduke/reference_data/archive/master.zip'
REFERENCE_DATA_DIR = os.path.join(ROTAMER_OUTPUT_DIR, 'reference_data-master', 'Top8000', 'Top8000_rotamer_pct_contour_grids')
FILENAMES = { 'ARG' : 'rota8000-arg.data',
              'ASN' : 'rota8000-asn.data',
              'ASP' : 'rota8000-asp.data',
              'CYS' : 'rota8000-cys.data',
              'GLN' : 'rota8000-gln.data',
              'GLU' : 'rota8000-glu.data',
              'HIS' : 'rota8000-his.data',
              'ILE' : 'rota8000-ile.data',
              'LEU' : 'rota8000-leu.data',
              'LYS' : 'rota8000-lys.data',
              'MET' : 'rota8000-met.data',
              'PHE' : 'rota8000-phetyr.data',
              'PRO' : 'rota8000-pro.data',
              'SER' : 'rota8000-ser.data',
              'THR' : 'rota8000-thr.data',
              'TRP' : 'rota8000-trp.data',
              'TYR' : 'rota8000-phetyr.data',
              'VAL' : 'rota8000-val.data' }
SERIALISED_LIB_PATH = os.path.join(ROTAMER_OUTPUT_DIR, 'library.pkl')
COMPRESSED_LIB_PATH = os.path.join(ROTAMER_OUTPUT_DIR, 'library.gz')
P_OUTLIER = 0.003
P_ALLOWED = 0.020


def download_reference_data():
    if not os.path.isdir(os.path.join(ROTAMER_OUTPUT_DIR, 'reference_data-master')):
        print('Downloading reference data...')
        response = requests.get(REFERENCE_DATA_URL, stream=True)
        with open(os.path.join(ROTAMER_OUTPUT_DIR, 'reference_data.zip'), 'wb') as outfile:
            for chunk in response.iter_content(chunk_size=1024):
                outfile.write(chunk)
        print('Unzipping reference data...')
        with zipfile.ZipFile(os.path.join(ROTAMER_OUTPUT_DIR, 'reference_data.zip'), 'r') as infile:
            infile.extractall(os.path.join(ROTAMER_OUTPUT_DIR))
        print('Deleting reference data ZIP...')
        os.remove(os.path.join(ROTAMER_OUTPUT_DIR, 'reference_data.zip'))
        print('Done.')


# These two functions are slow, hence why they were swapped out for NumPy bit-masking in the
# testing script and in the actual implementation. Don't use these for anything time critical.
# I only left these in here in case I want them for something else

# Compress a list of values into a list of integers of bit width n*l
def int_compress(xs, n, l):
    ys = [ ]
    xis = [ ]
    ps = [ 2**(l*(n-i-1)) for i in range(n) ]
    if len(xs)%n > 0:
        xs = xs + [ 0 ] * (n-len(xs)%n)
    for x in xs:
        xis.append(x)
        if len(xis) == n:
            y = 0
            for i, xi in enumerate(xis):
                y += xi * ps[i]
            ys.append(y)
            xis = [ ]
    return ys


# Decompress each value in a list into n integers of bit width l
def int_decompress(xs, n, l):
    ys = [ ]
    ps = [ 2**(l*(n-i-1)) for i in range(n) ]
    for x in xs:
        yis = [ ]
        for i in range(n):
            yi = x
            for j in range(i):
                yi -= yis[j] * ps[j]
            yi //= ps[i]
            yis.append(yi)
        ys += yis
    return ys


def generate_library():
    print('Generating rotamer library...')
    dim_offsets = { }
    dim_bin_ranges = { }
    dim_bin_widths = { }
    dim_num_options = { }
    compressed_classification_bytes = { }
    for code, filename in FILENAMES.items():
        print('*** ' + code)
        dim_offsets[code] = [ ]
        dim_bin_ranges[code] = [ ]
        dim_bin_widths[code] = [ ]
        dim_num_options[code] = [ ]
        # Load chi information and probabilities from file
        given_classifications = { }
        filepath = os.path.join(REFERENCE_DATA_DIR, filename)
        with open(filepath, 'r') as infile:
            for line in infile.readlines():
                if line[0] == '#':
                    if line[:5] == '#   x':
                        splitline = line.split(': ')[1].split(' ')
                        dim_min = float(splitline[0])
                        dim_max = float(splitline[1])
                        dim_size = dim_max - dim_min
                        dim_num_bins = int(splitline[2])
                        dim_bin_width = dim_size / float(dim_num_bins)
                        dim_bin_ranges[code].append((dim_min, dim_max))
                        dim_bin_widths[code].append(dim_bin_width)
                    continue
                splitline = [ float(x) for x in line.split(' ') ]
                chis = tuple(splitline[:-1])
                probability = float(splitline[-1])
                classification = 0
                if probability < P_OUTLIER:
                    classification = 1
                elif probability < P_ALLOWED:
                    classification = 2
                else:
                    classification = 3
                given_classifications[chis] = classification
        # Get first chi offsets
        for dimension, offset_chi in enumerate(given_classifications.keys()[0]):
            dim_range = dim_bin_ranges[code][dimension]
            dim_width = dim_bin_widths[code][dimension]
            new_offset_chi = offset_chi - dim_width * (offset_chi // dim_width)
            dim_offsets[code].append(new_offset_chi)
        # Get number of combinations and initialise the array
        dim_num_options[code] = [ int((r[1]-r[0])//w) for r, w in zip(dim_bin_ranges[code], dim_bin_widths[code]) ]
        total_num_combinations = np.product(dim_num_options[code])
        classifications = [ 0 ] * total_num_combinations
        # Insert known values into the array
        for chis, classification in given_classifications.items():
            index = 0
            for dimension, chi in enumerate(chis):
                dim_offset = dim_offsets[code][dimension]
                dim_width = dim_bin_widths[code][dimension]
                index += (chi - dim_offset) / dim_width * np.product(dim_num_options[code][dimension+1:])
            index = int(index)
            classifications[index] = classification
        # Generate compressed byte array
        classifications_compressed = int_compress(classifications, 4, 2)
        compressed_classification_bytes[code] = bytearray(classifications_compressed)
    # Write serialised data
    all_data = [ dim_offsets, dim_bin_ranges, dim_bin_widths, dim_num_options, compressed_classification_bytes ]
    with open(SERIALISED_LIB_PATH, 'wb') as outfile:
        pickle.dump(all_data, outfile)
    with gzip.open(COMPRESSED_LIB_PATH, 'wb') as outfile:
        pickle.dump(all_data, outfile)
    print('Done.')


if __name__ == '__main__':
    setup()
    download_reference_data()
    generate_library()
