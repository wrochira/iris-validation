import os

import matplotlib.pyplot as plt

from _defs import PDB_REDO_DATA_DIR, TIMING_OUTPUT_DIR
from common import setup, load_pdb_report_data, cleanup_all_pdb_redo_dirs


PDB_REPORT_DATA = { }


def draw_boxplots():
    data_original_mp = [ ]
    data_original_nomp = [ ]
    data_original_mp_per_res = [ ]
    data_original_nomp_per_res = [ ]
    with open(os.path.join(TIMING_OUTPUT_DIR, 'i2_results_0.csv'), 'r') as infile:
        infile.readline() # Skip header line
        for line in infile.readlines():
            splitline = line.strip().split(',')
            pdb_id = splitline[0]
            residue_count = int(PDB_REPORT_DATA[pdb_id]['residueCount'])
            num_errors = int(splitline[1])
            time_nomp_mean, time_nomp_sd, time_mp_mean, time_mp_sd = [ float(x) for x in splitline[2:] ]
            if num_errors == 0:
                data_original_mp.append(time_mp_mean)
                data_original_nomp.append(time_nomp_mean)
                data_original_mp_per_res.append(time_mp_mean / residue_count)
                data_original_nomp_per_res.append(time_nomp_mean / residue_count)

    data_new_mp = [ ]
    data_new_nomp = [ ]
    data_new_mp_per_res = [ ]
    data_new_nomp_per_res = [ ]
    with open(os.path.join(TIMING_OUTPUT_DIR, 'i2_results_1.csv'), 'r') as infile:
        infile.readline() # Skip header line
        for line in infile.readlines():
            splitline = line.strip().split(',')
            pdb_id = splitline[0]
            residue_count = int(PDB_REPORT_DATA[pdb_id]['residueCount'])
            num_errors = int(splitline[1])
            time_nomp_mean, time_nomp_sd, time_mp_mean, time_mp_sd = [ float(x) for x in splitline[2:] ]
            if num_errors == 0:
                data_new_mp.append(time_mp_mean)
                data_new_nomp.append(time_nomp_mean)
                data_new_mp_per_res.append(time_mp_mean / residue_count)
                data_new_nomp_per_res.append(time_nomp_mean / residue_count)

    plot_labels = [ 'Old Report', 'New Report\n(Iris and Molprobity)', 'New Report\n(Iris only)' ]
    plot_colors = [ 'moccasin', (0.8, 1.0, 0.8), (0.8, 1.0, 0.8) ]
    plot_values = [ [ data_original_mp, data_new_mp, data_new_nomp ],
                    [ data_original_mp_per_res, data_new_mp_per_res, data_new_nomp_per_res ] ]

    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 8))

    bplot1 = axes[0].boxplot(x=plot_values[0],
                             labels=plot_labels,
                             notch=True,
                             patch_artist=True,
                             showfliers=False,
                             widths=0.6)
    axes[0].set_ylim((0, None))
    axes[0].set_ylabel('Time (s)')
    axes[0].set_title('Run Time by Model')

    bplot2 = axes[1].boxplot(x=plot_values[1],
                             labels=plot_labels,
                             notch=True,
                             patch_artist=True,
                             showfliers=False,
                             widths=0.6)
    axes[1].set_ylim((0, None))
    axes[1].set_ylabel('Time (s)')
    axes[1].set_title('Run Time by Residue')

    for bplot in (bplot1, bplot2):
        for patch, color in zip(bplot['boxes'], plot_colors):
            patch.set_facecolor(color)

    plt.savefig(os.path.join(TIMING_OUTPUT_DIR, 'i2_boxplots.png'), dpi=300)


if __name__ == '__main__':
    setup()
    PDB_REPORT_DATA = load_pdb_report_data()
    draw_boxplots()
    cleanup_all_pdb_redo_dirs()
