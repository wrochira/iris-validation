"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

from math import pi
from shutil import rmtree
from os import mkdir, path

from iris_validation.metrics import get_percentile
from iris_validation.interface.charts import concentric, radar, grid
from iris_validation.utils import code_three_to_one, needleman_wunsch
from iris_validation import METRIC_NAMES, METRIC_POLARITIES, METRIC_SHORTNAMES, METRIC_DISPLAY_TYPES, REPORT_METRIC_IDS, REPORT_RESIDUE_VIEW

METRIC_NAMES = [ METRIC_NAMES[i] for i in REPORT_METRIC_IDS ]
METRIC_POLARITIES = [ METRIC_POLARITIES[i] for i in REPORT_METRIC_IDS ]
METRIC_SHORTNAMES = [ METRIC_SHORTNAMES[i] for i in REPORT_METRIC_IDS ]
METRIC_DISPLAY_TYPES = [ METRIC_DISPLAY_TYPES[i] for i in REPORT_METRIC_IDS ]

with open(path.join(path.dirname(path.realpath(__file__)), 'template', 'report.html'), 'r') as infile:
    HTML_TEMPLATE = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template', 'report_panel.html'), 'r') as infile:
    HTML_TEMPLATE_PANEL = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template', 'report_panel_bc.html'), 'r') as infile:
    HTML_TEMPLATE_PANEL_BC = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template', 'js', 'interaction.js'), 'r') as infile:
    JS_MINIFIED = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template', 'css', 'minimal.css'), 'r') as infile:
    CSS_MINIFIED = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template', 'css', 'compat.css'), 'r') as infile:
    CSS_COMPAT = infile.read()


def _align_chains(model_latest, model_previous):
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


def _generate_chain_sets(model_latest, model_previous):
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
        print('WARNING: at least one chain contains no amino acid residues. Ignoring chains: ' + ', '.join(sorted(bad_chain_ids)))
        for model in models:
            for chain_id in bad_chain_ids:
                model.remove_chain(chain_id)
    if 0 in [ model.chain_count for model in models ]:
        print('ERROR: no non-null chains in these models')
        return False
    return True


def _chart_data_from_models(model_latest, model_previous):
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
                metrics_continuous = [ residue.ramachandran_score, residue.rotamer_score, residue.avg_b_factor, residue.max_b_factor, residue.std_b_factor, residue.fit_score, residue.mainchain_fit_score, residue.sidechain_fit_score ]
                metrics_continuous = [ metrics_continuous[i] for i in REPORT_METRIC_IDS ]
                metrics_discrete = [ None for _ in REPORT_METRIC_IDS ]
                if 'Ramachandran Score' in METRIC_NAMES:
                    rama_index = METRIC_NAMES.index('Ramachandran Score')
                    metrics_discrete[rama_index] = None if residue.ramachandran_score is None else 2 if residue.ramachandran_favored else 1 if residue.ramachandran_allowed else 0
                if 'Rotamer Score' in METRIC_NAMES:
                    rota_index = METRIC_NAMES.index('Rotamer Score')
                    metrics_discrete[rota_index] = None if residue.rotamer_classification is None else 0 if residue.rotamer_classification in (0, 1) else residue.rotamer_classification-1
                marker = None
                if residue.molprobity_data is not None:
                    metrics_discrete[rota_index] = None if residue.rotamer_classification is None else 0 if residue.molprobity_data['rota_outlier'] else 2
                    marker = residue.molprobity_data['does_clash']
                metrics_percentiles = [ get_percentile(metric_id, value, model.resolution, normalise_polarity=True) for metric_id, value in zip(REPORT_METRIC_IDS, metrics_continuous) ]
                residue_chart_data = { 'continuous' : metrics_continuous,
                                       'discrete' : metrics_discrete,
                                       'marker' : marker,
                                       'percentiles' : metrics_percentiles,
                                       'code' : residue.code,
                                       'seqnum' : residue.sequence_number }
                chain_chart_data[pos_index][model_id] = residue_chart_data
        chart_data_by_cmr.append(chain_chart_data)
    return chart_data_by_cmr


def _concentric_charts_from_data(chart_data):
    molprobity_enabled = chart_data[0][0][-1]['marker'] is not None
    settings = { 'center_text_1' : 'Iris',
                 'center_text_2' : '',
                 'center_text_3' : '(Molprobity Enabled)' if molprobity_enabled else '',
                 'marker_label' : 'Clashes',
                 'ring_names' : METRIC_SHORTNAMES,
                 'ring_polarities' : METRIC_POLARITIES,
                 'ring_types' : METRIC_DISPLAY_TYPES,
                 'svg_hidden' : False,
                 'svg_id' : None }
    charts = [ ]
    for chain_id, chain_data in enumerate(chart_data):
        settings['center_text_2'] = 'Chain ' + chr(65+chain_id)
        settings['svg_hidden'] = False if chain_id == 0 else True
        settings['svg_id'] = 'iris-chart-' + str(chain_id)
        chart = concentric(chain_data, settings)
        charts.append(chart)
    return charts


def _radar_chart():
    settings = { 'axis_names' : METRIC_NAMES }
    return radar(settings)


def _grid_chart():
    settings = { 'box_1_label' : 'Ramachandran',
                 'box_2_label' : 'Rotamer',
                 'bar_label' : 'Percentiles',
                 'bar_1_label' : 'Avg. B-factor',
                 'bar_2_label' : 'Sidechain Fit' }
    return grid(settings)


def _generate_js_globals(chart_data):
    num_chains = len(chart_data)
    num_models = len(chart_data[0][0])
    num_metrics = len(chart_data[0][0][0])
    num_residues_by_chain = [ len(chain_data) for chain_data in chart_data ]

    js_string = '// Variables\n'
    js_string += 'let selectedModel = ' + str(num_models-1) + ';\n'
    js_string += 'let selectedChain = 0;\n'
    js_string += 'let selectedResidue = 0;\n'
    js_string += 'let isDragging = false;\n'
    js_string += 'let chartZoomed = false;\n'
    js_string += 'let modelMinMax = [ ];\n'
    js_string += 'let barOffsetY = 0;\n'
    js_string += 'let barMultiplierY = 0;\n'
    js_string += '// Constants\n'
    js_string += 'const numModels = ' + str(num_models) + ';\n'
    js_string += 'const numChains = ' + str(num_chains) + ';\n'
    js_string += 'const numMetrics = ' + str(num_metrics) + ';\n'
    js_string += 'const gapDegrees = ' + str(round(0.3*180/pi, 2)) + ';\n'
    js_string += 'const chainLengths = ' + str(num_residues_by_chain) + ';\n'

    mdim_js_strings = [ ]
    mdim_js_strings.append('const residueCodes = [ ')
    mdim_js_strings.append('const sequenceNumbers = [ ')
    mdim_js_strings.append('const absoluteMetrics = [ ')
    mdim_js_strings.append('const percentileMetrics = [ ')
    mdim_js_strings.append('const discreteMetrics = [ ')
    for model_id in range(num_models):
        for i in range(len(mdim_js_strings)):
            mdim_js_strings[i] += '[ '
        for chain_id in range(num_chains):
            for i in range(len(mdim_js_strings)):
                mdim_js_strings[i] += '[ '
            for dp_set in chart_data[chain_id]:
                residue = dp_set[model_id]
                if residue is None:
                    for i in range(len(mdim_js_strings)):
                        mdim_js_strings[i] += 'null, '
                    continue
                # Residue codes
                if residue['code'] is None:
                    mdim_js_strings[0] += 'null, '
                else:
                    mdim_js_strings[0] += '\'' + residue['code'] + '\'' + ', '
                # Sequence numbers
                if residue['seqnum'] is None:
                    mdim_js_strings[1] += 'null, '
                else:
                    mdim_js_strings[1] += str(residue['seqnum']) + ', '
                # Metric values - absolute
                if residue['continuous'] is None:
                    mdim_js_strings[2] += 'null, '
                else:
                    mdim_js_strings[2] += str([ None if x is None else round(x, 3) for x in residue['continuous'] ]).replace('None', 'null') + ', '
                # Metric values - percentiles
                if residue['percentiles'] is None:
                    mdim_js_strings[3] += 'null, '
                else:
                    mdim_js_strings[3] += str([ None if x is None else round(x, 3) for x in residue['percentiles'] ]).replace('None', 'null') + ', '
                # Metric values - discrete
                if residue['discrete'] is None:
                    mdim_js_strings[4] += 'null, '
                else:
                    mdim_js_strings[4] += str([ None if x is None else round(x, 3) for x in residue['discrete'] ]).replace('None', 'null') + ', '
            for i in range(len(mdim_js_strings)):
                mdim_js_strings[i] = mdim_js_strings[i][:-2] + ' ], '
        for i in range(len(mdim_js_strings)):
            mdim_js_strings[i] = mdim_js_strings[i][:-2] + ' ], '
    for i in range(len(mdim_js_strings)):
        mdim_js_strings[i] = mdim_js_strings[i][:-2] + ' ];\n'
    js_string += ''.join(mdim_js_strings)

    return js_string


def build_report(model_latest, model_previous, output_dir, mode=''):
    chart_data = _chart_data_from_models(model_latest, model_previous)

    concentric_charts = _concentric_charts_from_data(chart_data)
    if REPORT_RESIDUE_VIEW.lower() == 'grid':
        residue_chart = _grid_chart()
    elif REPORT_RESIDUE_VIEW.lower() == 'radar':
        residue_chart = _radar_chart()
    else:
        print('WARNING: Unrecognised report residue view mode specified in _defs. Recognised modes are "Grid" and "Radar"')
        residue_chart = ''

    ip1_html = ''
    for chain_id in range(len(chart_data)):
        colourDefinition = ' color: rgb(150,150,150);' if chain_id == 0 else ''
        ip1_html += '            <div style="float: left; margin-right: 20px;">\n'
        ip1_html += '              <a id="chain-button-' + str(chain_id) + '" onclick="setChain(' + str(chain_id) + ');" style="text-decoration: none;' + colourDefinition + '"><font size="4">Chain ' + chr(65+chain_id) + '</font></a>\n'
        ip1_html += '            </div>\n'
    ip2_html = ''.join([ '              ' + concentric_chart + '\n' for concentric_chart in concentric_charts ])[:-1] # [:-1] to remove trailing newline
    ip3_html = '              ' + residue_chart

    report_html = HTML_TEMPLATE_PANEL if mode == 'panel' else HTML_TEMPLATE_PANEL_BC if mode == 'panel-bc' else HTML_TEMPLATE
    report_html = report_html \
                    .replace('INJECTION_POINT_1', ip1_html) \
                    .replace('INJECTION_POINT_2', ip2_html) \
                    .replace('INJECTION_POINT_3', ip3_html)

    if path.isdir(output_dir):
        rmtree(output_dir)
    mkdir(output_dir)
    mkdir(path.join(output_dir, 'js'))
    mkdir(path.join(output_dir, 'css'))

    js_interaction = JS_MINIFIED
    js_globals = _generate_js_globals(chart_data)
    # The browser used in CCP4-i2 doesn't support the "let" statement
    if mode == 'panel-bc':
        js_interaction = js_interaction.replace('let ', 'var ')
        js_globals = js_globals.replace('let ', 'var ')

    with open(path.join(output_dir, 'report.html'), 'w') as outfile:
        outfile.write(report_html)
    with open(path.join(output_dir, 'js', 'interaction.js'), 'w') as outfile:
        outfile.write(js_interaction)
    with open(path.join(output_dir, 'js', 'globals.js'), 'w') as outfile:
        outfile.write(js_globals)
    with open(path.join(output_dir, 'css', 'minimal.css'), 'w') as outfile:
        outfile.write(CSS_MINIFIED)
    with open(path.join(output_dir, 'css', 'compat.css'), 'w') as outfile:
        outfile.write(CSS_COMPAT)
