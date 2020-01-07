from math import pi, cos, sin, ceil

import svgwrite
import numpy as np


COLOURS = { 'BLACK' : svgwrite.rgb(0, 0, 0),
			'L_GREY' : svgwrite.rgb(150, 150, 150),
			'D_GREY' : svgwrite.rgb(50, 50, 50),
			'WHITE' : svgwrite.rgb(255, 255, 255),
			'RED' : svgwrite.rgb(200, 50, 50),
			'GREEN' : svgwrite.rgb(50, 200, 50),
			'BLUE' : svgwrite.rgb(50, 50, 200),
			'YELLOW' : svgwrite.rgb(250, 250, 50),
			'ORANGE' : svgwrite.rgb(250, 200, 50),
			'MAGENTA' : svgwrite.rgb(200, 50, 200),
			'CYAN' : svgwrite.rgb(50, 200, 200) }

COLOURS_SEQUENCE = [ COLOURS['BLACK'], COLOURS['RED'], COLOURS['GREEN'], COLOURS['BLUE'], COLOURS['CYAN'], COLOURS['MAGENTA'], COLOURS['ORANGE'] ]
#COLOURS_SEQUENCE = [ COLOURS['BLACK'], COLOURS['BLUE'] ] * 10
#COLOURS_SEQUENCE = [ COLOURS['BLACK'] ] * 10

def _coords_from_angle(centre, angle, point_radius, adj=(0, 0), gap=0):
	xc, yc = centre
	x1 = xc + point_radius * sin(angle + gap/2) + adj[0]
	y1 = yc - point_radius * cos(angle + gap/2) + adj[1]
	return (round(x1, 2), round(y1, 2))


def concentric(protein_group, chain_id, sidepanel_height, canvas=(1000, 1000), gap=0.3):
	# Shortnames
	metric_names = protein_group.metric_names
	metric_shortnames = protein_group.metric_shortnames
	metric_stats_rep = protein_group.metric_stats_rep
	metric_stats_abs = protein_group.metric_stats_abs
	metric_quartiles_abs = protein_group.metric_quartiles_abs
	metric_sqdiff_minmax_rep = protein_group.metric_sqdiff_minmax_rep
	metric_sqdiff_minmax_abs = protein_group.metric_sqdiff_minmax_abs
	latest_chain = protein_group.proteins[-1].chains[chain_id]

	# Useful calculations
	sizes = (min(canvas)//100, min(canvas)//60)
	centre = (canvas[0]//2, canvas[1]//2)
	full_radius = min(canvas)//2 - 10
	division_size = full_radius // (len(metric_names)+2)

	# Initialise drawing
	dwg = svgwrite.Drawing(profile='full')

	# Set viewbox attribute for scaling
	dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in canvas ])
	dwg.attribs['width'] = '95%'
	dwg.attribs['height'] = '95%'

	# Set HTML attributes for interaction
	dwg.attribs['id'] = 'concentric-chart-' + str(chain_id)
	if chain_id > 0:
		dwg.attribs['style'] = 'display: none;'

	# For each protein, loop through each metric and draw the circle in fragments
	for protein_id, protein in enumerate(protein_group.proteins):
		chain = protein.chains[chain_id]
		angle_delta = (2*pi - gap) / len(chain.residues)
		for i in range(len(metric_names)):
			metric_mean, metric_sdev = metric_stats_rep[i]
			metric_min, metric_max = metric_sqdiff_minmax_rep[i]
			colour = COLOURS_SEQUENCE[i]
			metric_fragments = [ ]
			current_fragment = [ ]
			for j, residue in enumerate(chain.residues):
				metric_value = residue.metrics[i]['represented_value']
				if metric_value is None:
					if len(current_fragment) > 0:
						residue_point = _coords_from_angle(centre, angle_delta*(j+0.5), (i+2)*division_size, gap=gap)
						current_fragment.append(residue_point)
						metric_fragments.append(current_fragment)
						current_fragment = [ ]
						current_fragment.append(residue_point)
					continue
				metric_squared_diff = ((metric_value - metric_mean)/metric_sdev)**2
				polarity = 1 if metric_value >= metric_mean else -1
				polarised_squared_diff = polarity * metric_squared_diff
				point_radius = (i+2)*division_size + polarised_squared_diff/(metric_max-metric_min)*division_size*0.9
				residue_point = _coords_from_angle(centre, angle_delta*(j+0.5), point_radius, gap=gap)
				current_fragment.append(residue_point)
			if len(current_fragment) > 0:
				metric_fragments.append(current_fragment)
				current_fragment = [ ]
			#metric_fragments[-1].append(metric_fragments[0][0])
			if protein_id == len(protein_group.proteins)-1:
				for fragment in metric_fragments:
					dwg.add(dwg.polyline(points=fragment,
										fill_opacity=0,
										stroke=colour,
										stroke_width=2,
										stroke_opacity=1))
			else:
				for fragment in metric_fragments:
					dwg.add(dwg.polyline(points=fragment,
										fill=COLOURS['L_GREY'],
										fill_opacity=0,
										stroke=COLOURS['L_GREY'],
										stroke_width=2,
										stroke_opacity=1.0,
										stroke_dasharray=4))

	# For the main protein, draw markers for missing data and reset the fragment list
	for i in range(len(metric_names)):
		colour = COLOURS_SEQUENCE[i]
		for j, residue in enumerate(latest_chain.residues):
			metric_value = residue.metrics[i]['represented_value']
			if metric_value is None:
				residue_point = _coords_from_angle(centre, angle_delta*(j+0.5), (i+2)*division_size, gap=gap)
				dwg.add(dwg.circle(r=sizes[0]/1.8,
								   center=residue_point,
								   fill=colour,
								   fill_opacity=1,
								   stroke=COLOURS['WHITE'],
								   stroke_width=2,
								   stroke_opacity=1))

	# Chain markings around the circle
	dwg.add(dwg.circle(r=full_radius-2*sizes[1],
					   center=centre,
					   fill_opacity=0,
					   stroke=COLOURS['BLACK'],
					   stroke_width=1,
					   stroke_opacity=0.25))
	dwg.add(dwg.circle(r=full_radius-1*sizes[1],
					   center=centre,
					   fill_opacity=0,
					   stroke=COLOURS['BLACK'],
					   stroke_width=1,
					   stroke_opacity=0.25))
	for i in range(latest_chain.length+1):
		dwg.add(dwg.line(_coords_from_angle(centre, angle_delta*i, full_radius-2*sizes[1], gap=gap),
						 _coords_from_angle(centre, angle_delta*i, full_radius-1*sizes[1], gap=gap),
						 stroke=COLOURS['BLACK'],
						 stroke_width=1,
						 stroke_opacity=0.35))

	# Coloured chain markings around the circle based on each residue's metrics
	for i, residue in enumerate(latest_chain.residues):
		quartile_scores = [ ]
		for j, metric in enumerate(residue.metrics):
			metric_value = metric['absolute_value']
			metric_min, q1, q2, q3, metric_max = metric_quartiles_abs[j]
			if metric['optimisation'] == 'maximise':
				score = 0 if metric_value < q1 else 1 if metric_value < q2 else 2 if metric_value < q3 else 3
			if metric['optimisation'] == 'minimise':
				score = 3 if metric_value < q1 else 2 if metric_value < q2 else 1 if metric_value < q3 else 0
			quartile_scores.append(score)
		avg_score = int(ceil(sum(quartile_scores)/len(quartile_scores)))
		avg_colour = [ COLOURS['RED'], COLOURS['YELLOW'], COLOURS['GREEN'], COLOURS['GREEN'] ][avg_score]
		dwg.add(dwg.polygon(points=[ _coords_from_angle(centre, angle_delta*i, full_radius-2*sizes[1], gap=gap),
									 _coords_from_angle(centre, angle_delta*i, full_radius-1*sizes[1], gap=gap),
									 _coords_from_angle(centre, angle_delta*(i+1), full_radius-1*sizes[1], gap=gap),
									 _coords_from_angle(centre, angle_delta*(i+1), full_radius-2*sizes[1], gap=gap) ],
							stroke_opacity=0,
							fill=avg_colour,
							fill_opacity=0.2,
							onclick='setRadarChart(' + str(chain_id) + ',' + str(i) + ');'))

	# Draw white-out triangle for excess outer circles
	dwg.add(dwg.polygon(points=[ centre,
								 _coords_from_angle(centre, -gap/2, full_radius),
								 _coords_from_angle(centre, +gap/2, full_radius) ],
						stroke_opacity=0.25,
						fill=COLOURS['WHITE'],
						fill_opacity=1.0))
	for i in (0, latest_chain.length):
		dwg.add(dwg.line(centre,
						 _coords_from_angle(centre, angle_delta*i, full_radius-2*sizes[1], gap=gap),
						 stroke=COLOURS['BLACK'],
						 stroke_width=1,
						 stroke_opacity=0.25,
						 stroke_dasharray=0))

	# Draw axis labels
	for i, shortname in enumerate(metric_shortnames):
		dwg.add(dwg.polyline([ _coords_from_angle(centre, (gap/25)*(j-(20-1)/2), (i+2)*division_size) for j in range(20) ],
						 stroke=COLOURS_SEQUENCE[i],
						 stroke_width=1,
						 stroke_opacity=1,
						 fill_opacity=0,
						 stroke_dasharray=0))

		dwg.add(dwg.text(text=shortname,
						 insert=_coords_from_angle(centre, angle_delta*0, (i+2)*division_size+sizes[1]),
						 font_size=sizes[1],
						 font_family='Arial',
						 text_anchor='middle',
						 alignment_baseline='central'))

	# Label chain numbers around the circle
	for i in range(latest_chain.length):
		dwg.add(dwg.text(text=str(i+1),
						 insert=_coords_from_angle(centre, angle_delta*(i+0.5), full_radius, gap=gap),
						 font_size=sizes[1],
						 opacity=0,#1 if i in (0, sidepanel_height-1) else 0,
						 font_family='Arial',
						 text_anchor='middle',
						 alignment_baseline='central',
						 id='residue-label-' + str(latest_chain.chain_id) + '-' + str(i)))

	# Scrolling marker
	pass

	# Scrolling segment and overlay segments for interaction
	'''
	seglen = min(sidepanel_height, latest_chain.length)
	points = [ centre ] + [ _coords_from_angle(centre, i*angle_delta, full_radius-1*sizes[1], gap=gap) for i in range(seglen+1) ]
	dwg.add(dwg.polygon(points=points,
						fill=COLOURS['L_GREY'],
						fill_opacity=0.10,
						stroke=COLOURS['L_GREY'],
						stroke_opacity=1,
						id='segment' + str(latest_chain.chain_id)))
	for i in range(latest_chain.length-1):
		dwg.add(dwg.polygon(points=[ centre, _coords_from_angle(centre, angle_delta*i, full_radius, gap=gap), _coords_from_angle(centre, angle_delta*(i+1), full_radius, gap=gap) ],
							fill=COLOURS['L_GREY'],
							fill_opacity=0,
							stroke=COLOURS['D_GREY'],
							stroke_opacity=0,
							id='seg' + str(i),
							onmousedown='handleSegment(1, ' + str(latest_chain.chain_id) + ', ' + str(i) + ', ' + str(round(gap*180/pi, 2)) +');',
							onmouseover='handleSegment(2, ' + str(latest_chain.chain_id) + ', ' + str(i) + ', ' + str(round(gap*180/pi, 2)) +');',
							onmouseup='handleSegment(3, ' + str(latest_chain.chain_id) + ', ' + str(i) + ', ' + str(round(gap*180/pi, 2)) +');'))
	'''

	return dwg.tostring()


def concentric_whole_protein(protein, sidepanel_height, canvas=(1000, 1000)):
	# Shortnames
	chain_intervals = [ sum([ chain.length for chain in protein.chains[:i] ]) for i in range(len(protein.chains)) ]
	all_residues = protein.all_residues
	metric_names = protein.metric_names
	metric_stats_rep = protein.metric_stats_rep
	metric_sqdiff_minmax_rep = protein.metric_sqdiff_minmax_rep

	# Useful calculations
	sizes = (min(canvas)//100, min(canvas)//60)
	angle_delta = (2*pi) / len(all_residues)
	centre = (canvas[0]//2, canvas[1]//2)
	full_radius = min(canvas)//2 - 10
	division_size = full_radius // (len(metric_names)+1)

	# Initialise drawing
	dwg = svgwrite.Drawing(profile='full')

	# Set viewbox attribute for scaling
	dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in canvas ])
	dwg.attribs['width'] = '95%'
	dwg.attribs['height'] = '95%'

	# Set HTML attributes for interaction
	dwg.attribs['id'] = 'concentric-chart-whole-protein'
	dwg.attribs['style'] = 'display: none;'

	# Loop through each metric and draw the circle in fragments
	for i in range(len(metric_names)):
		metric_mean, metric_sdev = metric_stats_rep[i]
		metric_min, metric_max = metric_sqdiff_minmax_rep[i]
		colour = COLOURS_SEQUENCE[i]
		metric_fragments = [ ]
		current_fragment = [ ]
		for j, residue in enumerate(all_residues):
			metric_value = residue.metrics[i]['represented_value']
			if metric_value is None:
				if len(current_fragment) > 0:
					residue_point = _coords_from_angle(centre, angle_delta*(j+0.5), (i+1)*division_size)
					current_fragment.append(residue_point)
					metric_fragments.append(current_fragment)
					current_fragment = [ ]
					current_fragment.append(residue_point)
				continue
			metric_squared_diff = ((metric_value - metric_mean)/metric_sdev)**2
			polarity = 1 if metric_value >= metric_mean else -1
			polarised_squared_diff = polarity * metric_squared_diff
			point_radius = (i+1)*division_size + polarised_squared_diff/(metric_max-metric_min)*division_size*0.9
			residue_point = _coords_from_angle(centre, angle_delta*(j+0.5), point_radius)
			current_fragment.append(residue_point)
		if len(current_fragment) > 0:
			metric_fragments.append(current_fragment)
			current_fragment = [ ]
		metric_fragments[-1].append(metric_fragments[0][0])
		for fragment in metric_fragments:
			dwg.add(dwg.polyline(points=fragment,
								 fill=COLOURS['WHITE'],
								 fill_opacity=0,
								 stroke=colour,
								 stroke_width=2,
								 stroke_opacity=1))

	# Draw markers for missing data and reset the fragment list
	for i in range(len(metric_names)):
		colour = COLOURS_SEQUENCE[i]
		for j, residue in enumerate(all_residues):
			metric_value = residue.metrics[i]['represented_value']
			if metric_value is None:
				residue_point = _coords_from_angle(centre, angle_delta*(j+0.5), (i+1)*division_size)
				dwg.add(dwg.circle(r=sizes[0]/1.8,
								   center=residue_point,
								   fill=colour,
								   fill_opacity=1,
								   stroke=COLOURS['WHITE'],
								   stroke_width=2,
								   stroke_opacity=1))

	# Dividing lines to illustrate chain breaks
	for i in chain_intervals:
		dwg.add(dwg.line(centre,
						 _coords_from_angle(centre, angle_delta*i, (len(metric_names)+0.5)*division_size),
						 stroke=COLOURS['BLACK'],
						 stroke_width=2,
						 stroke_opacity=0.50,
						 stroke_dasharray=4))

	return dwg.tostring()


def scrolling_bar(protein, chain_id, sidepanel_height, canvas=(200, 1000)):
	# Shortnames
	chain = protein.chains[chain_id]
	metric_names = protein.metric_names
	metric_quartiles = protein.metric_quartiles_abs

	# Useful calculations
	bar_height = canvas[1]//sidepanel_height

	# Initialise drawing
	dwg = svgwrite.Drawing(profile='full')

	# Set viewbox attribute for scaling
	dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in canvas ])
	dwg.attribs['width'] = '95%'
	dwg.attribs['height'] = '95%'

	# Set HTML attributes for interaction
	dwg.attribs['id'] = 'bar-chart-' + str(chain.chain_id)
	if chain.chain_id > 0:
		dwg.attribs['style'] = 'display: none;'

	# Loop through each metric and draw each bar
	for i, metric_name in enumerate(metric_names):
		metric_min, q1, q2, q3, metric_max = metric_quartiles[i]
		for j, residue in enumerate(chain.residues):
			metric_value = residue.metrics[i]['absolute_value']
			optimisation = residue.metrics[i]['optimisation']
			if metric_value is None:
				continue
			if optimisation == 'maximise':
				normalised_metric_value = (metric_value-metric_min)/(metric_max-metric_min)
				bar_colour = COLOURS['RED'] if metric_value < q1 else COLOURS['ORANGE'] if metric_value < q2 else COLOURS['YELLOW'] if metric_value < q3 else COLOURS['GREEN']
			elif optimisation == 'minimise':
				normalised_metric_value = (metric_max-metric_value)/(metric_max-metric_min)
				bar_colour = COLOURS['GREEN'] if metric_value < q1 else COLOURS['YELLOW'] if metric_value < q2 else COLOURS['ORANGE'] if metric_value < q3 else COLOURS['RED']
			dwg.add(dwg.rect(insert=(0, j*bar_height),
							 size=(canvas[0]*normalised_metric_value, bar_height),
							 fill=bar_colour,
							 fill_opacity=0.2,
							 stroke=COLOURS['BLACK'],
							 stroke_width=1,
							 stroke_opacity=1,
							 onclick='setRadarChart(' + str(chain_id) + ',' + str(j) + ');'))

	return dwg.tostring()


def radar(protein, canvas=(600, 500)):
	# Shortnames
	metric_names = protein.metric_names

	# Useful calculations
	centre = (canvas[0]//2, canvas[1]//2)
	sizes = (min(canvas)//50, min(canvas)//40)
	angle_delta = 2*pi / len(metric_names)
	axis_radius = min(canvas)//2 - min(canvas)//20

	# Initialise drawing
	dwg = svgwrite.Drawing(profile='full')

	# Set viewbox attribute for scaling
	dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in canvas ])
	dwg.attribs['width'] = '100%'
	dwg.attribs['height'] = '100%'

	# Set HTML attributes for interaction
	dwg.attribs['id'] = 'radar-chart'

	# Axes
	for i in range(len(metric_names)):
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
	for i in range(len(metric_names)):
		l_width, l_height = sizes[1]/2 * len(metric_names[i]), sizes[1]
		vert_adj = l_height if pi/2 < angle_delta*i < 3*pi/2 else -l_height if 3*pi/2 < angle_delta*i or angle_delta*i >= 0 else 0
		horiz_adj = l_width/2 if 0 < angle_delta*i < pi else -l_width/2 if pi < angle_delta*i < 2*pi else 0
		dwg.add(dwg.text(text=metric_names[i],
						 insert=_coords_from_angle(centre, angle_delta*i, axis_radius, adj=(horiz_adj, vert_adj)),
						 font_size=sizes[1],
						 font_family='Arial',
						 text_anchor='middle',
						 alignment_baseline='central'))
	for i in range(10):
		label = str(10*(i+1))
		dwg.add(dwg.rect(insert=(centre[0] - sizes[0], centre[1] - (i+1)*axis_radius//10 - sizes[0]//2),
						 size=(sizes[0]*2, sizes[0]*1.5),
						 fill=COLOURS['WHITE']))
		dwg.add(dwg.text(text=label,
						 insert=(centre[0] - sizes[0]//1.5 - sizes[0]//3*(len(label)-2), centre[1] - (i+1)*axis_radius//10 + sizes[0]//2),
						 font_size=sizes[0],
						 font_family='Arial',
						 fill=COLOURS['D_GREY']))

	# Colour axes based on protein-wide metric percentiles
	'''
	chains_of_residues = [ chain.residues for chain in protein.chains ]
	all_residues = [ residue for chain in chains_of_residues for residue in chain ]
	protein_percentiles = _population_percentiles(metric_names, all_residues)
	for i, metric in enumerate(metric_names):
		mmin, mq1, mq2, mq3, mmax = [ x*axis_radius/100 for x in protein_percentiles[i] ]
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
	'''

	# Dots, labels, and connecting polygon (initialised hidden)
	dwg.add(dwg.polygon(points=[ _coords_from_angle(centre, angle_delta*i, axis_radius/2) for i in range(len(metric_names)) ],
						fill=COLOURS['WHITE'],
						fill_opacity=0,
						stroke=COLOURS['BLACK'],
						stroke_width=1,
						stroke_opacity=0,
						id='connecting-polygon'))
	for i in range(len(metric_names)):
		dwg.add(dwg.circle(r=sizes[0],
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
	floating_label.add(dwg.circle(r=sizes[1]*1.5,
								  center=(0, 0),
								  fill=COLOURS['WHITE'],
								  fill_opacity=0.8))
	floating_label.add(dwg.text('',
								insert=(0, 0),
								font_size=sizes[1],
								font_family='Arial',
								fill=COLOURS['BLACK'],
								fill_opacity=1,
								text_anchor='middle',
								alignment_baseline='central',
								id='floating-label-text'))
	dwg.add(floating_label)

	return dwg.tostring()
