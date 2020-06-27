"""
Copyright 2020 William Rochira at the University of York

Produces Iris-compatible data from Clipper-Tools MetricsModel objects
"""

import itertools

try:
    import clipper
    CLIPPER_MODE = 0
except ImportError:
    try:
        from clipper_python import _clipper as clipper
        CLIPPER_MODE = 1
    except ImportError:
        raise Exception('failed to import Clipper-Python.')

from clipper_tools.utils import code_three_to_one, needleman_wunsch
from clipper_tools.metrics import generate_metrics_model
from clipper_tools.metrics import get_percentile


def _align_chains(model_latest, model_previous=None):
    if model_previous is None:
        return
    latest_chain_ids = [ chain.chain_id for chain in model_latest ]
    previous_chain_ids = [ chain.chain_id for chain in model_previous ]
    extra_chains = set(latest_chain_ids) - set(previous_chain_ids)
    if len(extra_chains) > 0:
        raise NotImplementedError('Looks like you\'re finally going to have to write this code.')
    lost_chains = set(previous_chain_ids) - set(latest_chain_ids)
    if len(lost_chains) > 0:
        for lost_chain_id in lost_chains:
            model.remove_chain(lost_chain_id)
        print('WARNING: the chain count for the latest model is lower than that of previous models. These lost chains will not be represented in the validation report.')

def _clear_non_aa_residues(*models):
    for model in models:
       for chain in model:
            lost_residues = [ residue for residue in chain if not residue.is_aa ]
            for residue in lost_residues:
                chain.remove_residue(residue)

def _generate_chain_sets(model_latest, model_previous=None):
    chain_sets = { }
    for chain in model_latest:
        chain_sets[chain.chain_id] = [ chain ]
    if model_previous is not None:
        for chain in model_previous:
            if chain.chain_id in chain_sets:
                chain_sets[chain.chain_id] = [ chain, chain_sets[chain.chain_id][0] ]
    return chain_sets

def _get_chain_seqnum_minmax(chain_sets):
    minmax_by_chain = { }
    for chain_id, chain_set in chain_sets.items():
        minimum, maximum = None, None
        for chain in chain_set:
            for residue in chain:
                if minimum is None or residue.sequence_number < minimum:
                    minimum = residue.sequence_number
                if maximum is None or residue.sequence_number > maximum:
                    maximum = residue.sequence_number
        minmax_by_chain[chain_id] = (minimum, maximum)
    return minmax_by_chain

def _get_residue_alignments(chain_sets):
    # TODO: swap out Needleman-Wunsch for some MSA implementation if multiple-model alignment makes it into the release
    alignment_pair_by_chain = { }
    for chain_id, chain_set in chain_sets.items():
        sequences = [ code_three_to_one([ residue.code for residue in chain ]) for chain in chain_set ]
        if len(sequences) == 1:
            alignment_pair_by_chain[chain_id] = (sequences[0], )
            continue
        alignment_pair = needleman_wunsch(sequences[-2], sequences[-1])
        alignment_pair_by_chain[chain_id] = alignment_pair
    return alignment_pair_by_chain

def _verify_chain_lengths(*models):
    bad_chain_ids = set()
    for model in models:
        for chain in model:
            if chain.length == 0:
                bad_chain_ids.add(chain.chain_id)
    if len(bad_chain_ids) > 0:
        print('WARNING: at least one chain contains no amino acid residues. Removing chains: ' + ', '.join(bad_chain_ids))
        for model in models:
            for chain_id in bad_chain_ids:
                model.remove_chain(chain_id)
    if 0 in [ model.chain_count for model in models ]:
        print('ERROR: no non-null chains in these models')
        return False
    return True

def chart_data_from_models(model_latest, model_previous=None):
    models = [ model_latest, model_previous ]
    if model_previous is None:
        print('WARNING: previous model not supplied, chart will not support ghosting')
        models = [ model_latest ]

    _align_chains(*models)
    _clear_non_aa_residues(*models)
    cl_success = _verify_chain_lengths(*models)
    if not cl_success:
        return None
    chain_sets = _generate_chain_sets(*models)
    chain_seqnum_minmax = _get_chain_seqnum_minmax(chain_sets)
    alignment_pair_by_chain = _get_residue_alignments(chain_sets)

    chart_data_by_cmr = [ ] # CMR = chain, model, residue
    for chain_id, chain_set in chain_sets.items():
        aligned_length = len(alignment_pair_by_chain[chain_id][0])
        chain_chart_data = [ [ None for _ in range(len(models)) ] for _ in range(aligned_length) ]
        for model_id, model in enumerate(models):
            alignment = alignment_pair_by_chain[chain_id][model_id]
            residue_id = -1
            for pos_index, alignment_char in enumerate(alignment):
                if alignment_char == '-':
                    continue
                residue_id += 1
                residue = chain_set[model_id].residues[residue_id]
                metrics_continuous = [ residue.ramachandran_score, residue.rotamer_score, residue.avg_b_factor, residue.max_b_factor, residue.mainchain_fit_score, residue.sidechain_fit_score ]
                metrics_discrete = [ None, None, None, None, None, None ]
                metrics_discrete[0] = None if residue.ramachandran_score is None else 2 if residue.ramachandran_favored else 1 if residue.ramachandran_allowed else 0
                metrics_discrete[1] = None if residue.rotamer_classification is None else 0 if residue.rotamer_classification in (0, 1) else residue.rotamer_classification-1
                marker = None
                if residue.molprobity_data is not None:
                    metrics_discrete[1] = None if residue.rotamer_classification is None else 0 if residue.molprobity_data['rota_outlier'] else 2
                    marker = residue.molprobity_data['does_clash']
                metrics_percentiles = [ get_percentile(metric_id, value, model.resolution, normalise_polarity=True) for metric_id, value in enumerate(metrics_continuous) ]
                residue_chart_data = { 'continuous' : metrics_continuous,
                                       'discrete' : metrics_discrete,
                                       'marker' : marker,
                                       'percentiles' : metrics_percentiles,
                                       'code' : residue.code,
                                       'seqnum' : residue.sequence_number }
                chain_chart_data[pos_index][model_id] = residue_chart_data
        chart_data_by_cmr.append(chain_chart_data)
    return chart_data_by_cmr
