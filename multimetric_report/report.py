from os import mkdir, path

from shutil import rmtree

from . import charts, METRICS


with open(path.join(path.dirname(path.realpath(__file__)), 'template/report.html'), 'r') as infile:
	HTML_TEMPLATE = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/js/interaction.js'), 'r') as infile:
	JS_MINIFIED = infile.read()
with open(path.join(path.dirname(path.realpath(__file__)), 'template/css/minimal.css'), 'r') as infile:
	CSS_MINIFIED = infile.read()


def _generate_js_globals(protein):
	js_string = 'let lastClicked;\nlet chartZoomed = false;\n'
	js_string += 'const protein = [ '
	for chain in protein.chains:
		js_string += '[ '
		for residue in chain.residues:
			js_string += '{ \'aa\' : \'' + residue.aa + '\', \'metrics\' : ' + str(list(residue.metrics.values())) + ' }, '
		js_string = js_string[:-2] # To cut off the trailing ', '
		js_string += ' ], '
	js_string = js_string[:-2] # To cut off the trailing ', '
	js_string += ' ];\n'
	js_string += 'const metrics = [ '
	for metric in METRICS:
		js_string += '\'' + metric + '\'' + ', '
	js_string = js_string[:-2] # To cut off the trailing ', '
	js_string += ' ]'
	return js_string


def generate(protein, sidepanel_height, output_dir):
	concentric_charts = [ charts.concentric(chain, sidepanel_height) for chain in protein.chains ]
	radial_charts = [ charts.radial(chain) for chain in protein.chains ]
	radar_chart = charts.radar(protein)

	ip1_html = ''
	for i in range(len(protein.chains)):
		colourDefinition = ' color: rgb(150,150,150);' if i == 0 else ''
		ip1_html += '\t\t\t\t\t\t<div style="float: left; margin-right: 20px;">\n'
		ip1_html += '\t\t\t\t\t\t\t<a href="#" id="chain-button-' + str(i) + '" onclick="setChain(' + str(i) + ');" style="text-decoration: none;' + colourDefinition + '"><h6>Chain ' + str(i+1) + '</h6></a>\n'
		ip1_html += '\t\t\t\t\t\t</div>\n'
	#ip2_html = ''.join([ '\t\t\t\t\t\t\t' + radial_chart + '\n' for radial_chart in radial_charts ])
	ip2_html = ''.join([ '\t\t\t\t\t\t\t' + concentric_chart + '\n' for concentric_chart in concentric_charts ])
	ip3_html = '\t\t\t\t\t\t\t' + radar_chart + '\n'

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
		outfile.write(_generate_js_globals(protein))
