"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

import os
import sys
from math import exp, log, sqrt

try:
    import clipper
    CLIPPER_MODE = 0
except ImportError:
    try:
        from clipper_python import _clipper as clipper
        CLIPPER_MODE = 1
    except ImportError:
        raise ImportError('failed to import Clipper-Python')

from iris_validation.utils import ATOMIC_NUMBERS, MC_ATOM_NAMES, norm_cdf


class ReflectionsHandler(object):
    def __init__(self, f_reflections=None, xmap=None, header_mode=0):
        self.f_reflections = f_reflections
        self.xmap = xmap
        self.header_mode = header_mode

        self.hkl = clipper.HKL_info()
        self.f_phi = clipper.HKL_data_F_phi_float(self.hkl)

        self.headers_used = None
        self.grid = None
        self.spacegroup = None
        self.cell = None
        self.resolution = None
        self.resolution_limit = None

        if f_reflections is None:
            if xmap is None:
                print('ERROR: either a reflections file path or an xmap object must be passed as an argument')
                raise Exception('ArgumentError')
            try:
                self.grid = xmap.grid
            except AttributeError:
                self.grid = None
            self.spacegroup = xmap.spacegroup
            self.cell = xmap.cell
        else:
            extension = f_reflections.split('.')[-1].lower()
            if extension != 'mtz':
                if extension == 'cif':
                    print('ERROR: mmCIF format is not currently supported for reflections data.')
                else:
                    print('ERROR: reflections file has unrecognised extension:' + extension)
                raise Exception('ExtensionError')
            self._load_hkl_data()
            self._generate_xmap()
            self._calculate_map_stats()

    def _load_hkl_data(self):
        mtzin = clipper.CCP4MTZfile()
        try:
            mtzin.open_read(self.f_reflections)
            mtzin.import_hkl_info(self.hkl)
            header_hierarchy = [ ]
            # Need a better way to choose the right headers, but I'm not familiar enough with reflections data to know what labels are common
            """
            labels_and_types = [ tuple(str(line).strip().split(' ')) for line in mtzin.column_labels() ]
            labels_and_types_by_prefix = { }
            for column_label, column_type in labels_and_types:
                label_prefix = '/'.join(column_label.split('/')[:-1])
                label_suffix = column_label.split('/')[-1]
                if label_prefix not in labels_and_types_by_prefix:
                    labels_and_types_by_prefix[label_prefix] = [ ]
                labels_and_types_by_prefix[label_prefix].append((column_label, column_type))
            candidate_pair = None
            for prefix in labels_and_types_by_prefix.keys():
                column_labels, column_types = zip(*labels_and_types_by_prefix[prefix])
                F_labels, P_labels, Q_labels = [ ], [ ], [ ]
                for column_label, column_type in labels_and_types_by_prefix[prefix]:
                    if column_type == 'F':
                        F_labels.append(column_label)
                    if column_type == 'P':
                        P_labels.append(column_label)
                    if column_type == 'Q':
                        Q_labels.append(column_label)
            if candidate_pair is None:
                print('ERROR: reflections file does not contain the required columns')
            """
            if self.header_mode in (0, 2):
                header_hierarchy += [ ('F', 'PHI'),
                                      ('FC', 'PHIC'),
                                      ('FC_ALL', 'PHIC_ALL') ]
            if self.header_mode in (0, 1):
                header_hierarchy += [ ('F', 'SIGF'),
                                      ('FP', 'SIGFP'),
                                      ('FP_ALL', 'SIGFP_ALL') ]
            if CLIPPER_MODE == 0:
                mtz_labels_and_types = [ tuple(str(line).strip().split(' ')) for line in mtzin.column_labels() ]
            elif CLIPPER_MODE == 1:
                mtz_labels_and_types = [ tuple(str(line).strip().split(' ')) for line in mtzin.column_labels ]
            mtz_column_labels, mtz_column_types = zip(*mtz_labels_and_types)
            mtz_column_label_suffixes = set([ label.split('/')[-1] for label in mtz_column_labels ])
            import_status = -1
            for set_id, suffixes in enumerate(header_hierarchy):
                if len(set(suffixes) - mtz_column_label_suffixes) == 0:
                    try:
                        mtzin.import_hkl_data(self.f_phi, '/*/*/[' + ','.join(suffixes) + ']')
                        import_status = set_id 
                        break
                    except Exception as e:
                        print('ERROR: failed to import HKL data from reflections file')
                        raise e
            if import_status == -1:
                print('ERROR: reflections file does not contain the required columns')
                raise Exception('ColumnError')
            else:
                self.headers_used = ', '.join(header_hierarchy[import_status])
            mtzin.close_read()
        except Exception as e:
            print('ERROR: failed to load data from reflections file:', self.f_reflections)
            raise(e)

        if CLIPPER_MODE == 0:
            spacegroup = self.hkl.spacegroup()
            cell = self.hkl.cell()
            resolution = self.hkl.resolution()
        elif CLIPPER_MODE == 1:
            spacegroup = self.hkl.spacegroup
            cell = self.hkl.cell
            resolution = self.hkl.resolution

        self.spacegroup = spacegroup
        self.cell = cell
        self.resolution = resolution
        self.resolution_limit = self.resolution.limit()

    def _generate_xmap(self):
        self.grid = clipper.Grid_sampling(self.spacegroup, self.cell, self.resolution)
        self.xmap = clipper.Xmap_float(self.spacegroup, self.cell, self.grid)
        try:
            self.xmap.fft_from(self.f_phi)
        except AttributeError:
            self.xmap.fft_from_float(self.f_phi)

    def _calculate_map_stats(self):
        map_stats = clipper.Map_stats(self.xmap)
        self.map_mean = map_stats.mean()
        self.map_std = map_stats.std_dev()

    def get_density_at_point(self, xyz):
        if CLIPPER_MODE == 0:
            cell = self.xmap.cell()
            grid = self.xmap.grid_sampling()
        elif CLIPPER_MODE == 1:
            cell = self.xmap.cell
            grid = self.xmap.grid_sampling
        co = clipper.Coord_orth(*xyz)
        co_cf_cg = co.coord_frac(cell).coord_grid(grid)
        density = self.xmap.get_data(co_cf_cg)
        return density

    def get_density_at_atom(self, mmol_atom):
        if CLIPPER_MODE == 0:
            co = mmol_atom.coord_orth()
            xyz = (co.x(), co.y(), co.z())
        elif CLIPPER_MODE == 1:
            xyz = mmol_atom.coord
        return self.get_density_at_point(xyz)

    def get_density_scores_at_residue(self, metrics_residue):
        all_atom_scores, mainchain_atom_scores, sidechain_atom_scores = [ ], [ ], [ ]
        for atom_id, atom in enumerate(metrics_residue.minimol_residue):
            is_mainchain = str(atom.name()).strip() in MC_ATOM_NAMES
            if CLIPPER_MODE == 0:
                element = str(atom.element()).strip()
            elif CLIPPER_MODE == 1:
                element = atom.element.strip()
            atomic_number = ATOMIC_NUMBERS[element]
            density = self.get_density_at_atom(atom)
            #density_norm = density / atomic_number
            #atom_score = log(norm_cdf((density_norm - self.map_mean) / self.map_std))
            atom_score = -log(norm_cdf((density - self.map_mean) / self.map_std) / atomic_number)
            all_atom_scores.append(atom_score)
            if is_mainchain:
                mainchain_atom_scores.append(atom_score)
            else:
                sidechain_atom_scores.append(atom_score)
        all_score, mainchain_score, sidechain_score = None, None, None
        if len(all_atom_scores) > 0:
            all_atom_score_avg = sum(all_atom_scores) / len(all_atom_scores)
            all_score = all_atom_score_avg * metrics_residue.avg_b_factor
        if metrics_residue.is_aa:
            if len(mainchain_atom_scores) > 0:
                mainchain_atom_score_avg = sum(mainchain_atom_scores) / len(mainchain_atom_scores)
                mainchain_score = mainchain_atom_score_avg * metrics_residue.mc_b_factor
            if len(sidechain_atom_scores) > 0:
                sidechain_atom_score_avg = sum(sidechain_atom_scores) / len(sidechain_atom_scores)
                sidechain_score = sidechain_atom_score_avg * metrics_residue.sc_b_factor
        return all_score, mainchain_score, sidechain_score
