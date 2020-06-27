"""
Copyright 2020 William Rochira at the University of York
"""

from math import pi
from os import mkdir, path
from shutil import rmtree

from .charts import iris, radar, residue


METRIC_NAMES = ('Ramachandran Score', 'Rotamer Score', 'Avg B-factor', 'Max B-factor', 'Mainchain Fit', 'Sidechain Fit')
METRIC_POLARITIES = (+1, -1, -1, -1, -1, -1)
METRIC_SHORTNAMES = ('Rama', 'Rota', 'Avg B', 'Max B', 'MC Fit', 'SC Fit')
METRIC_DISPLAY_TYPES = ('D', 'D', 'C', 'C', 'C', 'C')


with open(path.join(path.dirname(path.realpath(__file__)), 'template/report.html'), 'r') as infile:
    HTML_TEMPLATE = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/report_minimal.html'), 'r') as infile:
    HTML_TEMPLATE_MINIMAL = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/report_pane.html'), 'r') as infile:
    HTML_TEMPLATE_PANE = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/js/interaction.js'), 'r') as infile:
    JS_MINIFIED = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/css/minimal.css'), 'r') as infile:
    CSS_MINIFIED = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/css/compat.css'), 'r') as infile:
    CSS_COMPAT = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/css/pane.css'), 'r') as infile:
    CSS_PANE = infile.read()


def _load_css_minimal():
    return CSS_MINIFIED

def _load_css_compat():
    return CSS_COMPAT

def _load_js_interaction():
    return JS_MINIFIED

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

def _iris_charts_from_data(chart_data):
    molprobity_enabled = chart_data[0][0][-1]['marker'] is not None
    iris_settings = { 'ring_names' : METRIC_SHORTNAMES,
                      'ring_types' : METRIC_DISPLAY_TYPES,
                      'ring_polarities' : METRIC_POLARITIES,
                      'center_text_1' : 'Iris',
                      'center_text_2' : '',
                      'center_text_3' : '(Molprobity enabled)' if molprobity_enabled else '',
                      'marker_label' : 'Clashes',
                      'svg_id' : None,
                      'svg_hidden' : False }
    iris_charts = [ ]
    for chain_id, chain_data in enumerate(chart_data):
        iris_settings['center_text_2'] = 'Chain ' + chr(65+chain_id)
        iris_settings['svg_hidden'] = False if chain_id == 0 else True
        iris_settings['svg_id'] = 'iris-chart-' + str(chain_id)
        iris_chart = iris(chain_data, iris_settings)
        iris_charts.append(iris_chart)
    return iris_charts

def _residue_chart():
    return residue()

def report_from_data(chart_data, output_dir, mode='', js_old=False):
    iris_charts = _iris_charts_from_data(chart_data)
    residue_chart = _residue_chart()

    ip1_html = ''
    for chain_id in range(len(chart_data)):
        colourDefinition = ' color: rgb(150,150,150);' if chain_id == 0 else ''
        ip1_html += '            <div style="float: left; margin-right: 20px;">\n'
        #ip1_html += '              <a href="#" id="chain-button-' + str(chain_id) + '" onclick="setChain(' + str(chain_id) + ');" style="text-decoration: none;' + colourDefinition + '"><font size="3">Chain ' + chr(65+chain_id) + '</font></a>\n'
        ip1_html += '              <a id="chain-button-' + str(chain_id) + '" onclick="setChain(' + str(chain_id) + ');" style="text-decoration: none;' + colourDefinition + '"><font size="3">Chain ' + chr(65+chain_id) + '</font></a>\n'
        ip1_html += '            </div>\n'
    ip2_html = ''.join([ '              ' + iris_chart + '\n' for iris_chart in iris_charts ])[:-1] # [:-1] to remove trailing newline
    #ip3_html = '              ' + radar_chart
    ip3_html = '              ' + residue_chart

    report_html = HTML_TEMPLATE_MINIMAL if mode == 'minimal' else HTML_TEMPLATE_PANE if mode == 'pane' else HTML_TEMPLATE
    report_html = report_html \
                    .replace('INJECTION_POINT_1', ip1_html) \
                    .replace('INJECTION_POINT_2', ip2_html) \
                    .replace('INJECTION_POINT_3', ip3_html)

    if path.isdir(output_dir):
        rmtree(output_dir)
    mkdir(output_dir)
    mkdir(path.join(output_dir, 'js'))
    mkdir(path.join(output_dir, 'css'))

    js_interaction = _load_js_interaction()
    js_globals = _generate_js_globals(chart_data)
    if js_old:
        js_interaction = js_interaction.replace('let ', 'var ')
        js_globals = js_globals.replace('let ', 'var ')

    with open(path.join(output_dir, 'report.html'), 'w') as outfile:
        outfile.write(report_html)
    with open(path.join(output_dir, 'js', 'interaction.js'), 'w') as outfile:
        outfile.write(js_interaction)
    with open(path.join(output_dir, 'js', 'globals.js'), 'w') as outfile:
        outfile.write(js_globals)
    if mode == 'pane':
        with open(path.join(output_dir, 'css', 'pane.css'), 'w') as outfile:
            outfile.write(CSS_PANE)
    else:
        with open(path.join(output_dir, 'css', 'minimal.css'), 'w') as outfile:
            outfile.write(CSS_MINIFIED)
        with open(path.join(output_dir, 'css', 'compat.css'), 'w') as outfile:
            outfile.write(CSS_COMPAT)
