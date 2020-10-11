"""
Copyright 2020 William Rochira at York Structural Biology Laboratory

- Simplifies calculation of per-residue metrics from a Clipper-Python MiniMol type
- Introduces types: MetricsModel, MetricsChain, MetricsResidue;
- MetricsModel must be initialised with a Clipper-Python MiniMol model type
- MetricsChain must be initialised with a Clipper-Python MiniMol polymer (chain) type
- MetricsResidue must be initialised with a Clipper-Python MiniMol monomer (residue) type, but also accepts arguments that contextualise it within the chain, allowing for further
  metric calculations, such as mainchain torsion angles and Ramachandran Plot probability
- When initialised by a MetricsChain instance, each MetricsResidue object will have all the context arguments specified

TODO:
- Integrate Molprobity analysis
"""

import math

from iris_validation import utils
from iris_validation import clipper
from iris_validation import SC_INCOMPLETE_STRING
from iris_validation.metrics.rotamer import get_cv_score
from iris_validation.metrics.percentiles import get_percentile
from iris_validation.metrics.reflections import ReflectionsHandler


#RAMA_THRESHOLDS = (0.01, 0.0005) # Clipper default
RAMA_THRESHOLDS = (0.02, 0.002) # Concordant with Coot


class MetricsModel(object):
    def __init__(self, mmol_model, reflections_handler=None, multithreaded=True):
        self._index = -1
        self.minimol_model = mmol_model
        self.reflections_handler = reflections_handler
        self.multithreaded = multithreaded

        minimol_chains = list(mmol_model.model())
        self.chain_count = len(minimol_chains)
        self.resolution = reflections_handler.resolution_limit if reflections_handler is not None else None

        # Multiprocessing is infeasible because of SWIG bindings
        # Multithreading tests revealed insignificant improvements in run time
        """
        from multiprocessing.pool import ThreadPool
        pool = ThreadPool(self.chain_count)
        self.chains = pool.map(_mmc_to_mmc, minimol_chains)
        """
        self.chains = [ MetricsChain(c, self) for c in mmol_model.model() ]

    def __iter__(self):
        return self

    # Python 3: def __next__(self):
    def next(self):
        if self._index < len(self.chains)-1:
            self._index += 1
            return self.chains[self._index]
        self._index = -1
        raise StopIteration

    def add_molprobity_data(self, clash_residues=[ ], rotamer_outliers=[ ]):
        #clash_residues = set(clash_residues)
        #rotamer_outliers = set(rotamer_outliers)
        for chain in self.chains:
            for residue in chain.residues:
                rid = (chain.chain_id, residue.sequence_number)
                clash = rid in clash_residues
                rota = rid in rotamer_outliers
                rmd = { 'does_clash' : clash,
                        'rota_outlier' : rota }
                residue.molprobity_data = rmd

    def get_chain(self, chain_id):
        return next(chain for chain in self.chains if chain.chain_id == chain_id)

    def remove_chain(self, chain_id):
        matching_chains = [ chain for chain in self.chains if chain.chain_id == chain_id ]
        if len(matching_chains) == 0:
            print('Error removing chain, no chains matching that ID were found.')
        else:
            for chain in matching_chains:
                self.chains.remove(chain)
                self.chain_count -= 1

    def b_factor_lists(self):
        all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs = [ [ ] for _ in range(8) ]
        for chain in self.chains:
            all_bfs_c, aa_bfs_c, mc_bfs_c, sc_bfs_c, non_aa_bfs_c, water_bfs_c, ligand_bfs_c, ion_bfs_c = chain.b_factor_lists()
            for model_li, chain_li in zip((all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs),
                                          (all_bfs_c, aa_bfs_c, mc_bfs_c, sc_bfs_c, non_aa_bfs_c, water_bfs_c, ligand_bfs_c, ion_bfs_c)):
                model_li += chain_li
        return all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs


class MetricsChain(object):
    def __init__(self, mmol_chain, parent=None):
        self.parent = parent
        self.minimol_chain = mmol_chain
        self._index = -1
        self.residues = [ ]
        self.length = len(mmol_chain)
        self.chain_id = str(mmol_chain.id().trim())
        for i, mmres in enumerate(mmol_chain):
            previous = mmol_chain[i-1] if i>0 else None
            next = mmol_chain[i+1] if i<len(mmol_chain)-1 else None
            residue = MetricsResidue(mmres, i, previous, next, self)
            self.residues.append(residue)
        for i, residue in enumerate(self.residues):
            if (i > 0 and i < len(self.residues)-1) and \
               (self.residues[i-1].is_aa and residue.is_aa and self.residues[i+1].is_aa) and \
               (self.residues[i-1].sequence_number+1 == residue.sequence_number == self.residues[i+1].sequence_number-1):
                residue.is_consecutive_aa = True
            else:
                residue.is_consecutive_aa = False

    def __iter__(self):
        return self

    # Python 3: def __next__(self):
    def next(self):
        if self._index < self.length-1:
            self._index += 1
            return self.residues[self._index]
        self._index = -1
        raise StopIteration

    def get_residue(self, sequence_number):
        return next(residue for residue in self.residues if residue.sequence_number == sequence_number)

    def remove_residue(self, residue):
        if residue in self.residues:
            self.residues.remove(residue)
            self.length -= 1
        else:
            print('Error removing residue, no matching residue was found.')

    def b_factor_lists(self):
        all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs = [ [ ] for _ in range(8) ]
        for residue in self.residues:
            all_bfs.append(residue.avg_b_factor)
            if residue.is_aa:
                aa_bfs.append(residue.avg_b_factor)
                mc_bfs.append(residue.mc_b_factor)
                sc_bfs.append(residue.sc_b_factor)
            else:
                non_aa_bfs.append(residue.avg_b_factor)
                if residue.is_water:
                    water_bfs.append(residue.avg_b_factor)
                # Silly assumption, just followed to be consistent with original i2 validation tool:
                elif len(residue.atoms) > 1:
                    ligand_bfs.append(residue.avg_b_factor)
                else:
                    ion_bfs.append(residue.avg_b_factor)
        return all_bfs, aa_bfs, mc_bfs, sc_bfs, non_aa_bfs, water_bfs, ligand_bfs, ion_bfs


class MetricsResidue(object):
    def __init__(self, mmol_residue, index_in_chain=None, previous=None, next=None, parent=None):
        self.parent = parent
        self.minimol_residue = mmol_residue
        self.initialised_with_context = index_in_chain is not None
        self.index_in_chain = index_in_chain
        self.previous = previous
        self.next = next
        self.atoms = [ atom for atom in mmol_residue ]
        self.sequence_number = mmol_residue.seqnum()
        self.code = mmol_residue.type().trim()
        self.code_type = utils.code_type(mmol_residue)
        self.backbone_atoms = utils.get_backbone_atoms(mmol_residue)
        self.backbone_atoms_are_correct = None not in self.backbone_atoms
        self.backbone_geometry_is_correct = utils.check_backbone_geometry(mmol_residue) if self.backbone_atoms_are_correct else None
        self.is_aa = utils.check_is_aa(mmol_residue)
        self.is_water = str(mmol_residue.type()).strip() == 'HOH'
        self.is_consecutive_aa = None
        self.phi = clipper.MMonomer.protein_ramachandran_phi(self.previous, mmol_residue) if self.previous else None
        self.psi = clipper.MMonomer.protein_ramachandran_psi(mmol_residue, self.next) if self.next else None
        if self.phi is not None and math.isnan(self.phi):
            self.phi = None
        if self.psi is not None and math.isnan(self.psi):
            self.psi = None
        self.chis = utils.calculate_chis(mmol_residue)
        self.is_sidechain_complete = SC_INCOMPLETE_STRING not in self.chis
        self.ramachandran_score = utils.calculate_ramachandran_score(mmol_residue, self.code, self.phi, self.psi)
        self.ramachandran_favored = RAMA_THRESHOLDS[0] <= self.ramachandran_score
        self.ramachandran_allowed = RAMA_THRESHOLDS[1] <= self.ramachandran_score < RAMA_THRESHOLDS[0]
        self.ramachandran_outlier = self.ramachandran_score < RAMA_THRESHOLDS[1]
        self.rotamer_score = utils.calculate_rotamer_score(mmol_residue, self.code, self.chis) if self.is_sidechain_complete else None
        self.rotamer_classification, self.rotamer_favored, self.rotamer_allowed = None, None, None
        if self.is_sidechain_complete:
            self.rotamer_classification = utils.get_rotamer_classification(mmol_residue, self.code, self.chis)
            self.rotamer_favored = True if self.rotamer_classification in (-1, 2) else False
            self.rotamer_allowed = False if self.rotamer_classification == 1 else False
            self.rotamer_outlier = True if self.rotamer_classification == 0 else False
        self.max_b_factor, self.avg_b_factor, self.std_b_factor, self.mc_b_factor, self.sc_b_factor = utils.analyse_b_factors(mmol_residue, self.is_aa, self.backbone_atoms)
        self.fit_score, self.mainchain_fit_score, self.sidechain_fit_score = None, None, None
        reflections_handler = parent.parent.reflections_handler
        if reflections_handler is not None:
            self.fit_score, self.mainchain_fit_score, self.sidechain_fit_score = reflections_handler.get_density_scores_at_residue(self)
        self.molprobity_data = None


def _mmc_to_mmc(minimol_chain):
    return MetricsChain(minimol_chain)


def generate_metrics_model(f_model=None, f_reflections=None, minimol=None, xmap=None, multithreaded=True):
    if f_model is None:
        if minimol is None:
            print('ERROR: either a model file path or a MiniMol object must be passed as an argument')
            return
    else:
        fpdb = clipper.MMDBfile()
        minimol = clipper.MiniMol()
        try:
            fpdb.read_file(f_model)
            fpdb.import_minimol(minimol)
        except Exception as e:
            print('ERROR: failed to import model file:')
            print(e)
            return
    if f_reflections is None:
        if xmap is None:
            print('WARNING: neither a reflections file path nor a map were specified; continuing in geometry-only mode')
            reflections_handler = None
        else:
            reflections_handler = ReflectionsHandler(xmap=xmap, minimol=minimol)
    else:
        reflections_handler = ReflectionsHandler(f_reflections, minimol=minimol)

    metrics_model = MetricsModel(minimol, reflections_handler, multithreaded)
    return metrics_model
