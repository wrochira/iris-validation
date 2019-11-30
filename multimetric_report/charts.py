from math import pi, cos, sin

import svgwrite
import numpy as np

from . import METRICS


COLOURS = { 'BLACK' : svgwrite.rgb(0, 0, 0),
			'L_GREY' : svgwrite.rgb(150, 150, 150),
			'D_GREY' : svgwrite.rgb(50, 50, 50),
			'WHITE' : svgwrite.rgb(255, 255, 255),
			'RED' : svgwrite.rgb(200, 50, 50),
			'GREEN' : svgwrite.rgb(50, 200, 50),
			'BLUE' : svgwrite.rgb(50, 50, 200),
			'YELLOW' : svgwrite.rgb(200, 200, 50),
			'MAGENTA' : svgwrite.rgb(200, 50, 200),
			'CYAN' : svgwrite.rgb(50, 200, 200) }


def _coords_from_angle(centre, angle, point_radius, adj=(0, 0), gap=0):
	xc, yc = centre
	x1 = xc + point_radius * sin(angle + gap/2) + adj[0]
	y1 = yc - point_radius * cos(angle + gap/2) + adj[1]
	return (round(x1, 2), round(y1, 2))


def _population_percentiles(residues):
	stats = { }
	for metric in METRICS:
		values = [ residue.metrics[metric] for residue in residues ]
		stats[metric] = [ round(np.percentile(values, q), 2) for q in (0, 25, 50, 75, 100) ]
	return stats


def concentric(chain, sidepanel_height, canvas=(1000, 1000)):
	# Useful calculations
	font_sizes = (min(canvas)//100, min(canvas)//60)
	angle_delta = (2*pi) / len(chain.residues)
	centre = (canvas[0]//2, canvas[1]//2)
	full_radius = min(canvas)//2 - 10
	division_size = full_radius // (len(METRICS)+1)

	# Initialise drawing
	dwg = svgwrite.Drawing(profile='full')

	# Set viewbox attribute for scaling
	dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in canvas ])
	dwg.attribs['width'] = '95%'
	dwg.attribs['height'] = '95%'

	# Set HTML attributes for interaction
	dwg.attribs['id'] = 'radial-chart-' + str(chain.chain_id)
	if chain.chain_id > 0:
		dwg.attribs['style'] = 'display: none;'

	# Baseline circles
	for i in range(len(METRICS)):
		dwg.add(dwg.circle(center=centre,
						   r=division_size*(i+1),
						   fill_opacity=0,
						   stroke=COLOURS['BLACK'],
						   stroke_width=1,
						   stroke_opacity=0.8))

	# Loop through each metric to analyse relative squared-deviation distributions
	for i, metric in enumerate(METRICS):
		all_metric_values = [ residue.metrics[metric] for residue in chain.residues ]
		mean_metric_value = sum(all_metric_values) / len(all_metric_values)
		#modal_metric_values = gaussian_kde(all_metric_values)

	# Loop through each metric and draw the circle
	for i, metric in enumerate(METRICS):
		metric_points = [ ]
		for j, residue in enumerate(chain.residues):
			metric_squared_diff = (residue.metrics[metric] - mean_metric_value)**2
			polarity = 1 if residue.metrics[metric] >= mean_metric_value else -1
			polarised_squared_diff = polarity * metric_squared_diff
			point_radius = (i+1)*division_size + max(-5000, min(5000, polarised_squared_diff))/5000 * (division_size)
			residue_point = _coords_from_angle(centre, angle_delta*j, point_radius)
			metric_points.append(residue_point)
		colour = [ COLOURS['BLACK'], COLOURS['RED'], COLOURS['GREEN'], COLOURS['BLUE'], COLOURS['CYAN'], COLOURS['YELLOW'], COLOURS['MAGENTA'] ][i]
		#metric_xys = np.asfortranarray(list(zip(*metric_points)))
		#smooth_points = bezier.Curve(metric_xys, degree=2)
		dwg.add(dwg.polygon(points=metric_points,
							fill=COLOURS['WHITE'],
							fill_opacity=0,
							stroke=colour,
							stroke_width=2,
							stroke_opacity=1))

	# Scrolling segment and overlay segments for interaction
	seglen = min(sidepanel_height, len(chain.residues))
	points = [ centre ] + [ _coords_from_angle(centre, (i-seglen/2)*angle_delta, full_radius) for i in range(seglen+1) ]
	dwg.add(dwg.polygon(points=points,
						fill=COLOURS['L_GREY'],
						fill_opacity=0.4,
						stroke=COLOURS['WHITE'],
						stroke_opacity=1,
						id='segment' + str(chain.chain_id)))
	for i in range(len(chain.residues)):
		dwg.add(dwg.polygon(points=[centre, _coords_from_angle(centre, angle_delta*i, full_radius), _coords_from_angle(centre, angle_delta*(i+1), full_radius)],
							fill=COLOURS['L_GREY'],
							fill_opacity=0,
							stroke=COLOURS['D_GREY'],
							stroke_opacity=0,
							id='seg' + str(i),
							onclick='handleSegment(1, ' + str(chain.chain_id) + ', ' + str(i) +');',
							onmouseover='handleSegment(2, ' + str(chain.chain_id) + ', ' + str(i) +');',
							onmouseout='handleSegment(3, ' + str(chain.chain_id) + ', ' + str(i) +');'))

	return dwg.tostring()


def radar(protein, canvas=(600, 500)):
	# Useful calculations
	centre = (canvas[0]//2, canvas[1]//2)
	font_sizes = (min(canvas)//50, min(canvas)//40)
	angle_delta = 2*pi / len(METRICS)
	axis_radius = min(canvas)//2 - min(canvas)//20

	# Initialise drawing
	dwg = svgwrite.Drawing(profile='full')

	# Set viewbox attribute for scaling
	dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in canvas ])
	dwg.attribs['width'] = '95%'
	dwg.attribs['height'] = '95%'

	# Set HTML attributes for interaction
	dwg.attribs['id'] = 'radar-chart'

	# Axes
	for i in range(len(METRICS)):
		dwg.add(dwg.line(centre,
						 _coords_from_angle(centre, angle_delta*i, axis_radius),
						 stroke=COLOURS['L_GREY'],
						 stroke_width=2))
		for j in range(10):
			dwg.add(dwg.line(_coords_from_angle(centre, angle_delta*i, (j+1)*axis_radius//10),
							 _coords_from_angle(centre, angle_delta*(i+1), (j+1)*axis_radius//10),
							 stroke=COLOURS['L_GREY'],
							 stroke_width=1))

	# Axis labels
	for i in range(len(METRICS)):
		l_width, l_height = font_sizes[1]/2 * len(METRICS[i]), font_sizes[1]
		vert_adj = l_height if pi/2 < angle_delta*i < 3*pi/2 else -l_height if 3*pi/2 < angle_delta*i or angle_delta*i >= 0 else 0
		horiz_adj = l_width/2 if 0 < angle_delta*i < pi else -l_width/2 if pi < angle_delta*i < 2*pi else 0
		dwg.add(dwg.text(text=METRICS[i],
						 insert=_coords_from_angle(centre, angle_delta*i, axis_radius, adj=(horiz_adj, vert_adj)),
						 font_size=font_sizes[1],
						 font_family='Arial',
						 text_anchor='middle',
						 alignment_baseline='central'))
	for i in range(10):
		label = str(10*(i+1))
		dwg.add(dwg.rect(insert=(centre[0] - font_sizes[0], centre[1] - (i+1)*axis_radius//10 - font_sizes[0]//2),
						 size=(font_sizes[0]*2, font_sizes[0]*1.5),
						 fill=COLOURS['WHITE']))
		dwg.add(dwg.text(text=label,
						 insert=(centre[0] - font_sizes[0]//1.5 - font_sizes[0]//3*(len(label)-2), centre[1] - (i+1)*axis_radius//10 + font_sizes[0]//2),
						 font_size=font_sizes[0],
						 font_family='Arial',
						 fill=COLOURS['D_GREY']))

	# Colour axes based on protein-wide metric percentiles
	chains_of_residues = [ chain.residues for chain in protein.chains ]
	all_residues = [ residue for chain in chains_of_residues for residue in chain ]
	protein_percentiles = _population_percentiles(all_residues)
	for i, metric in enumerate(METRICS):
		mmin, mq1, mq2, mq3, mmax = [ x*axis_radius/100 for x in protein_percentiles[metric] ]
		dwg.add(dwg.line(_coords_from_angle(centre, angle_delta*i, mmin),
						 _coords_from_angle(centre, angle_delta*i, mq1),
						 stroke=COLOURS['RED'],
						 stroke_width=4,
						 stroke_opacity=0.5))
		dwg.add(dwg.line(_coords_from_angle(centre, angle_delta*i, mq1),
						 _coords_from_angle(centre, angle_delta*i, mq3),
						 stroke=COLOURS['YELLOW'],
						 stroke_width=12,
						 stroke_opacity=0.5))
		dwg.add(dwg.line(_coords_from_angle(centre, angle_delta*i, mq3),
						 _coords_from_angle(centre, angle_delta*i, mmax),
						 stroke=COLOURS['GREEN'],
						 stroke_width=4,
						 stroke_opacity=0.5))

	# Dots, labels, and connecting polygon (initialised hidden)
	dwg.add(dwg.polygon(points=[ _coords_from_angle(centre, angle_delta*i, axis_radius/2) for i in range(len(METRICS)) ],
						fill=COLOURS['WHITE'],
						fill_opacity=0,
						stroke=COLOURS['BLACK'],
						stroke_width=1,
						stroke_opacity=0,
						id='connecting-polygon'))
	for i in range(len(METRICS)):
		dwg.add(dwg.circle(r=font_sizes[0],
						   #center=centre,
						   center=(0, 0),
						   fill=COLOURS['L_GREY'],
						   fill_opacity=0,
						   stroke=COLOURS['BLACK'],
						   stroke_width=1,
						   stroke_opacity=0,
						   id='metric-point-' + str(i),
						   onclick='handleRadar(1, ' + str(i) +');',
						   onmouseover='handleRadar(2, ' + str(i) +');',
						   onmouseout='handleRadar(3, ' + str(i) +');'))
	floating_label = dwg.g(id='floating-label-group',
						   opacity=0)
	floating_label.add(dwg.circle(r=font_sizes[1]*1.5,
								  center=(0, 0),
								  fill=COLOURS['WHITE'],
								  fill_opacity=0.8))
	floating_label.add(dwg.text('',
								insert=(0, 0),
								font_size=font_sizes[1],
								font_family='Arial',
								fill=COLOURS['BLACK'],
								fill_opacity=1,
								text_anchor='middle',
								alignment_baseline='central',
								id='floating-label-text'))
	dwg.add(floating_label)

	return dwg.tostring()


def radial(chain, canvas=(1000, 1000), gap=0.2):
	# Useful calculations
	middle_radius = min(canvas) // 5
	font_sizes = (min(canvas)//100, min(canvas)//60)
	angle_delta = (2*pi - gap) / len(chain.residues)
	centre = (canvas[0]//2, canvas[1]//2)
	division_size = (min(canvas) - 2 * middle_radius) // 20 - 1 # '-1' introduces 10px border
	full_radius = min(canvas)//2 - 10

	# Initialise drawing
	dwg = svgwrite.Drawing(profile='full')

	# Set viewbox attribute for scaling
	dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in canvas ])
	dwg.attribs['width'] = '95%'
	dwg.attribs['height'] = '95%'

	# Set HTML attributes for interaction
	dwg.attribs['id'] = 'radial-chart-' + str(chain.chain_id)
	if chain.chain_id > 0:
		dwg.attribs['style'] = 'display: none;'

	# Add integrated JS functions
	#dwg.add(dwg.script(content='function settext(doc_id, text) { document.getElementById(doc_id).innerHTML = text; }'))

	# Axes
	for i in range(11):
		dwg.add(dwg.circle(center=centre,
						   r=middle_radius+i*division_size+1,
						   fill_opacity=0,
						   stroke=COLOURS['BLACK'],
						   stroke_width=2,
						   stroke_opacity=0.8))

	# Radial stacked traffic lights
	for i, residue in enumerate(chain.residues):
		metric_values = sorted([ v for k, v in residue.metrics.items() ], reverse=True)
		for metric_value in metric_values:
			colour = COLOURS['RED'] if metric_value < 33 else COLOURS['YELLOW'] if metric_value < 67 else COLOURS['GREEN']
			metric_radius = middle_radius + (division_size * metric_value / 10)
			dwg.add(dwg.polygon(fill_opacity=0.5,
								points=[centre,
										_coords_from_angle(centre, angle_delta*i, metric_radius, gap=gap),
										_coords_from_angle(centre, angle_delta*(i+1), metric_radius, gap=gap)],
								fill=colour,
								stroke=COLOURS['WHITE'],
								stroke_width=1))
			# Line: dwg.add(dwg.line(centre, _coords_from_angle(centre, angle_delta*(i+0.5), metric_radius, gap=gap), stroke=colour, stroke_width=5))
			# Point: dwg.add(dwg.circle(r=2, center=_coords_from_angle(centre, angle_delta*(i+0.5), metric_radius, gap=gap), stroke=colour, stroke_width=5))

	# RMS line for all metrics
	for i, residue in enumerate(chain.residues):
		metric_values = sorted([ v for k, v in residue.metrics.items() ], reverse=True)
		# TODO: calculate RMS and draw thick line, maybe move this into the traffic light section

	# Centre circle
	dwg.add(dwg.circle(center=centre, r=middle_radius, fill=COLOURS['WHITE']))

	# 1-letter codes
	for i, residue in enumerate(chain.residues):
		dwg.add(dwg.polygon(points=[centre, _coords_from_angle(centre, angle_delta*i, middle_radius, gap=gap), _coords_from_angle(centre, angle_delta*(i+1), middle_radius, gap=gap)],
							fill=COLOURS['L_GREY'] if i%2 else COLOURS['D_GREY'],
							fill_opacity=0.2))
		dwg.add(dwg.text(residue.aa,
						 insert=_coords_from_angle(centre, angle_delta*(i+0.5), middle_radius-10, adj=(-font_sizes[0]//2, font_sizes[0]//2), gap=gap),
						 font_size=font_sizes[0],
						 font_family='Arial',
						 fill=COLOURS['BLACK']))
		dwg.add(dwg.text(i+1,
						 insert=_coords_from_angle(centre, angle_delta*(i+0.5), middle_radius-25, adj=(-font_sizes[0]//2, font_sizes[0]//2), gap=gap),
						 font_size=font_sizes[0],
						 font_family='Arial',
						 fill=COLOURS['BLACK']))
		'''dwg.add(dwg.line(_coords_from_angle(centre, angle_delta*(i+0.5), 0, gap=gap),
						 _coords_from_angle(centre, angle_delta*(i+0.5), middle_radius-40, gap=gap),
						 stroke=COLOURS['BLACK'],
						 stroke_width=1))'''

	# Axis labels
	dwg.add(dwg.polygon(fill_opacity=1,
						points=[centre,
								_coords_from_angle(centre, -gap/4, full_radius+10),
								_coords_from_angle(centre, gap/4, full_radius+10)],
						fill=COLOURS['WHITE'],
						stroke_opacity=0))
	for i in range(11):
		label = str(10*i)
		dwg.add(dwg.text(label,
						 insert=(centre[0] - font_sizes[1]//1.5 - font_sizes[1]//3*(len(label)-2), centre[1] - middle_radius - i*division_size + font_sizes[1]//3),
						 font_size=font_sizes[1],
						 font_family='Arial',
						 fill=COLOURS['D_GREY']))

	# Top segments for interaction
	for i in range(len(chain.residues)):
		dwg.add(dwg.polygon(points=[centre, _coords_from_angle(centre, angle_delta*i, full_radius, gap=gap), _coords_from_angle(centre, angle_delta*(i+1), full_radius, gap=gap)],
							fill=COLOURS['L_GREY'],
							fill_opacity=0,
							stroke=COLOURS['D_GREY'],
							stroke_opacity=0,
							id='seg' + str(i),
							onclick='handleChain(1, ' + str(chain.chain_id) + ', ' + str(i) +');',
							onmouseover='handleChain(2, ' + str(chain.chain_id) + ', ' + str(i) +');',
							onmouseout='handleChain(3, ' + str(chain.chain_id) + ', ' + str(i) +');'))

	return dwg.tostring()
