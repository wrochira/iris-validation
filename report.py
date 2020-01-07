from os import mkdir, path

from shutil import rmtree

from . import charts


with open(path.join(path.dirname(path.realpath(__file__)), 'template/report.html'), 'r') as infile:
	HTML_TEMPLATE = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/js/interaction.js'), 'r') as infile:
	JS_MINIFIED = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/css/minimal.css'), 'r') as infile:
	CSS_MINIFIED = infile.read()


def _generate_js_globals(protein, sidepanel_height):
	js_string = 'let lastClicked;\nlet chartZoomed = false;\nlet segmentDocked = true;\nlet sidepanelHeight = ' + str(sidepanel_height) + ';\n'
	js_string += 'const protein = [ '
	for chain in protein.chains:
		js_string += '[ '
		for residue in chain.residues:
			js_string += '{ \'code\' : \'' + residue.code + '\', \'metrics-absolute\' : ' + str([ metric['absolute_value'] for metric in residue.metrics ]).replace('None', 'null') + ', \'metrics-represented\' : ' + str([ metric['represented_value'] for metric in residue.metrics ]).replace('None', 'null') + ' }, '
		js_string = js_string[:-2] # To cut off the trailing ', '
		js_string += ' ], '
	js_string = js_string[:-2] # To cut off the trailing ', '
	js_string += ' ];\n'
	js_string += 'const totalResidues = ' + str(len(protein.all_residues)) + ';'
	return js_string


def generate(protein_group, sidepanel_height, output_dir):
	latest_protein = protein_group.proteins[-1]
	concentric_charts = [ charts.concentric(protein_group, i, sidepanel_height) for i in range(len(latest_protein.chains)) ] + [ charts.concentric_whole_protein(latest_protein, sidepanel_height) ]
	radar_chart = charts.radar(latest_protein)

	ip1_html = ''
	for i in range(len(latest_protein.chains)):
		colourDefinition = ' color: rgb(150,150,150);' if i == 0 else ''
		ip1_html += '\t\t\t\t\t\t<div style="float: left; margin-right: 20px;">\n'
		ip1_html += '\t\t\t\t\t\t\t<a href="#" id="chain-button-' + str(i) + '" onclick="setConcentric(' + str(i) + ');" style="text-decoration: none;' + colourDefinition + '"><h6>Chain ' + chr(65+i) + '</h6></a>\n'
		ip1_html += '\t\t\t\t\t\t</div>\n'
	ip1_html += '\t\t\t\t\t\t<div style="float: left; margin-right: 20px;">\n'
	ip1_html += '\t\t\t\t\t\t\t<a href="#" id="chain-button-whole-protein" onclick="setConcentric(\'whole-protein\');" style="text-decoration: none;"><h6>Whole Protein</h6></a>\n'
	ip1_html += '\t\t\t\t\t\t</div>\n'
	ip2_html = ''.join([ '\t\t\t\t\t\t\t\t' + concentric_chart + '\n' for concentric_chart in concentric_charts ])[:-1]
	ip3_html = '\t\t\t\t\t\t\t' + radar_chart

	report_html = HTML_TEMPLATE \
					.replace('INJECTION_POINT_1', ip1_html) \
					.replace('INJECTION_POINT_2', ip2_html) \
					.replace('INJECTION_POINT_3', ip3_html)

	if path.isdir(output_dir):
		rmtree(output_dir)
	mkdir(output_dir)
	mkdir(path.join(output_dir, 'js'))
	mkdir(path.join(output_dir, 'css'))

	with open(path.join(output_dir, 'report.html'), 'w') as outfile:
		outfile.write(report_html)
	with open(path.join(output_dir, 'css', 'minimal.css'), 'w') as outfile:
		outfile.write(CSS_MINIFIED)
	with open(path.join(output_dir, 'js', 'interaction.js'), 'w') as outfile:
		outfile.write(JS_MINIFIED)
	with open(path.join(output_dir, 'js', 'globals.js'), 'w') as outfile:
		outfile.write(_generate_js_globals(latest_protein, sidepanel_height))
