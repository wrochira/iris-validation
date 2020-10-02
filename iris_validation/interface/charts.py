"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

from math import pi, cos, sin

import svgwrite
from svgwrite.animate import Animate
from svgwrite.gradients import LinearGradient


COLORS = { 'BLACK' : 'rgb(0, 0, 0)',
           'VL_GREY' : 'rgb(200, 200, 200)',
           'L_GREY' : 'rgb(150, 150, 150)',
           'D_GREY' : 'rgb(50, 50, 50)',
           'WHITE' : 'rgb(255, 255, 255)',
           'RED' : 'rgb(200, 80, 80)',
           'L_RED' : 'rgb(220, 150, 150)',
           'VL_RED' : 'rgb(255, 235, 235)',
           'GREEN' : 'rgb(50, 200, 50)',
           'L_GREEN' : 'rgb(150, 220, 150)',
           'BLUE' : 'rgb(50, 50, 200)',
           'L_BLUE' : 'rgb(150, 150, 220)',
           'YELLOW' : 'rgb(250, 250, 50)',
           'L_YELLOW' : 'rgb(255, 255, 150)',
           'ORANGE' : 'rgb(250, 200, 50)',
           'L_ORANGE' : 'rgb(255, 220, 150)',
           'MAGENTA' : 'rgb(200, 50, 200)',
           'L_MAGENTA' : 'rgb(220, 150, 220)',
           'CYAN' : 'rgb(50, 200, 200)',
           'L_CYAN' : 'rgb(150, 220, 220)',
           'BAR_GREEN' : 'rgb(90, 237, 141)',
           'BAR_ORANGE' : 'rgb(247, 212, 134)',
           'BAR_RED' : 'rgb(240, 106, 111) '}

DEFAULT_SETTINGS_CONCENTRIC = { 'animation_length' : 250,
                                'apply_markers' : True,
                                'axis_min' : -0.7,
                                'axis_max' : 0.25,
                                'baseline_point_resolution' : 200,
                                'canvas_size' : (1000, 1000),
                                'center_text_1' : None,
                                'center_text_2' : None,
                                'center_text_3' : None,
                                'discrete_color_sequence' : ([ COLORS['RED'], COLORS['ORANGE'], COLORS['GREEN'] ]),
                                'discrete_ring_color' : COLORS['L_GREY'],
                                'top_gap': 0.3,
                                'interaction_segments_enabled' : True,
                                'marker_label' : '',
                                'marker_color' : COLORS['RED'],
                                'missing_data_shade_enabled' : True,
                                'missing_data_shade_color' : COLORS['VL_RED'],
                                'ring_color_sequence' : (COLORS['CYAN'], COLORS['MAGENTA'], COLORS['BLUE'], COLORS['ORANGE'], COLORS['CYAN'], COLORS['MAGENTA'], COLORS['BLUE'], COLORS['ORANGE']),
                                'ring_names' : None,
                                'ring_polarities' : None,
                                'ring_types' : None,
                                'segment_function_name' : 'handleSegment',
                                'segment_selector_enabled' : True,
                                'svg_id' : 'concentric',
                                'svg_hidden' : None }

DEFAULT_SETTINGS_RADAR = { 'canvas_size' : (600, 500),
                           'svg_id' : 'radar',
                           'svg_hidden' : None }

DEFAULT_SETTINGS_GRID = { 'box_1_label' : '',
                          'box_2_label' : '',
                          'bar_label' : '',
                          'bar_1_label' : '',
                          'bar_2_label' : '',
                          'canvas_size' : (400, 1000),
                          'svg_id' : 'grid',
                          'svg_hidden' : None }

_CFA_CACHE = { }


def _coords_from_angle(center, angle, point_radius, adj=(0, 0), gap=0):
    global _CFA_CACHE
    arg_string = str([ center, angle, point_radius, adj, gap ])
    if arg_string in _CFA_CACHE:
        coords = _CFA_CACHE[arg_string]
    else:
        xc, yc = center
        x1 = xc + point_radius * sin(angle + gap/float(2)) + adj[0]
        y1 = yc - point_radius * cos(angle + gap/float(2)) + adj[1]
        coords = (round(x1, 1), round(y1, 1))
        _CFA_CACHE[arg_string] = coords
    return coords


def concentric(datapoints, settings={ }):
    # Check arguments
    try:
        assert len(datapoints) > 0
        for dp_set in datapoints:
            for dp in dp_set:
                if dp is None:
                    continue
                if len(set(['continuous', 'discrete' ]) - set(dp.keys())) > 0:
                    print('ERROR: every datapoint must have at least "continuous" and "discrete" key values defined')
                    return ''
    except:
        print('ERROR: datapoints argument should be a two-dimensional iterable of dictionaries')
        return ''

    num_segments = len(datapoints)
    num_versions = len(datapoints[0])
    num_rings = len(datapoints[0][0]['continuous'])

    for setting_name, setting_default in DEFAULT_SETTINGS_CONCENTRIC.items():
        if setting_name not in settings.keys():
            settings[setting_name] = setting_default

    if settings['ring_names'] is None:
        settings['ring_names'] = [ '' for _ in range(num_rings) ]
    if len(settings['ring_names']) != num_rings:
        print('ERROR: specified number of ring names is inconsistent with the data')
        return ''
    if settings['ring_types'] is None:
        settings['ring_types'] = [ 'C' for _ in range(num_rings) ]
    if len(settings['ring_types']) != num_rings:
        print('ERROR: specified number of ring types is inconsistent with the data')
        return ''
    if settings['ring_polarities'] is None:
        settings['ring_polarities'] = [ +1 for _ in range(num_rings) ]
    if len(settings['ring_polarities']) != num_rings:
        print('ERROR: specified number of ring polarities is inconsistent with the data')
        return ''

    # Useful calculations
    gap = float(settings['top_gap'])
    sizes = (min(settings['canvas_size'])//100, min(settings['canvas_size'])//60)
    center = (settings['canvas_size'][0]//2, settings['canvas_size'][1]//2)
    full_radius = min(settings['canvas_size'])//2 - 10
    division_size = full_radius // (num_rings+2)
    angle_delta = (2*pi - gap) / float(num_segments)

    # Initialise drawing
    dwg = svgwrite.Drawing(profile='full')

    # Set viewbox attribute for scaling
    dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in settings['canvas_size'] ])
    dwg.attribs['width'] = '95%'
    dwg.attribs['height'] = '95%'

    # Set HTML attributes for interaction
    # If the user doesn't set an ID for the chart, then it shouldn't have an ID defined in the XML. But, the svg_id variable still needs to be
    # defined for later, to make sure there's something to prepend to the IDs of any child elements that *do* need IDs.
    if settings['svg_id'] is None:
        settings['svg_id'] = 'concentric'
    else:
        dwg.attribs['id'] = settings['svg_id']
    if settings['svg_hidden']:
        dwg.attribs['style'] = 'display: none;'

    # Draw axes
    for ring_id, ring_name in enumerate(settings['ring_names']):
        ring_color = settings['ring_color_sequence'][ring_id] if settings['ring_types'][ring_id] == 'C' else settings['discrete_ring_color']
        dwg.add(dwg.circle(r=(ring_id+2)*division_size,
                           center=center,
                           fill_opacity=0,
                           stroke=ring_color,
                           stroke_width=1,
                           stroke_opacity=1))
        dwg.add(dwg.polyline([ _coords_from_angle(center, (gap/float(25))*(i-(20-1)/float(2)), (ring_id+2)*division_size) for i in range(20) ],
                         stroke=ring_color,
                         stroke_width=3,
                         stroke_opacity=1,
                         fill_opacity=0))
        dwg.add(dwg.text(text=ring_name,
                         insert=_coords_from_angle(center, 0, (ring_id+2)*division_size+sizes[1]),
                         font_size=sizes[1],
                         font_family='Arial',
                         text_anchor='middle',
                         alignment_baseline='central'))

    # Get means for continuous metrics; get ranges for discrete metrics
    ring_averages = [ ]
    ring_categories = [ ]
    for ring_id in range(num_rings):
        ring_values_cont = [ ]
        ring_values_disc = set()
        for dp_set in datapoints:
            for dp in dp_set:
                if dp is None:
                    continue
                if dp['continuous'] is not None and dp['continuous'][ring_id] is not None:
                    dp['continuous'][ring_id] *= settings['ring_polarities'][ring_id]
                    ring_values_cont.append(dp['continuous'][ring_id])
                if dp['discrete'] is not None and dp['discrete'][ring_id] is not None:
                    ring_values_disc.add(dp['discrete'][ring_id])
        ring_avg = None
        if len(ring_values_cont) > 0:
            ring_avg = sum(ring_values_cont) / len(ring_values_cont)
        ring_averages.append(ring_avg)
        ring_ordered_set = tuple(sorted(ring_values_disc))
        ring_categories.append(ring_ordered_set)

    # Calculate magnitudes for every point and use them to find the min/max for each metric
    magnitudes_by_rsv = [ ] # RSV = ring, segment, version
    plot_magnitudes_by_rsv = [ ]
    for ring_id in range(num_rings):
        # Calculate average negative delta in the latest dataset
        latest_negative_deltas = [ ]
        for dp_set in datapoints:
            if dp_set[-1] is None or dp_set[-1]['continuous'][ring_id] is None:
                continue
            delta = dp_set[-1]['continuous'][ring_id] - ring_averages[ring_id]
            if delta < 0:
                latest_negative_deltas.append(delta)
        avg_negative_delta = 0
        if len(latest_negative_deltas) > 0:
            avg_negative_delta = sum(latest_negative_deltas) / len(latest_negative_deltas)
        # Subtract the average negative delta from all deltas to calculate 'magnitudes'
        ring_magnitudes = [ ]
        for dp_set in datapoints:
            magnitudes_set = [ ]
            for dp in dp_set:
                magnitude = None
                if dp is None or dp['continuous'] is None or dp['continuous'][ring_id] is None:
                    magnitudes_set.append(magnitude)
                    continue
                delta = dp['continuous'][ring_id] - ring_averages[ring_id]
                magnitude = delta - avg_negative_delta
                magnitudes_set.append(magnitude)
            ring_magnitudes.append(magnitudes_set)
        magnitudes_by_rsv.append(ring_magnitudes)
        # Calculate 'plot magnitudes'
        all_ring_magnitudes = [ magnitude for magnitudes_set in ring_magnitudes for magnitude in magnitudes_set if magnitude is not None ]
        ring_magnitude_minmax = (0, 0)
        if len(all_ring_magnitudes) > 0:
            ring_magnitude_minmax = (min(all_ring_magnitudes), max(all_ring_magnitudes))
        ring_plot_magnitudes = [ ]
        for dp_id, dp_set in enumerate(datapoints):
            plot_magnitudes_set = [ ]
            for version_id, dp in enumerate(dp_set):
                magnitude = ring_magnitudes[dp_id][version_id]
                if magnitude is None:
                    plot_magnitude = None
                else:
                    plot_magnitude = 0
                    if magnitude > 0 and ring_magnitude_minmax[1] != 0:
                        plot_magnitude = magnitude / ring_magnitude_minmax[1] * settings['axis_max']
                    elif magnitude < 0 and ring_magnitude_minmax[0] != 0:
                        plot_magnitude = magnitude / ring_magnitude_minmax[0] * settings['axis_min']
                plot_magnitudes_set.append(plot_magnitude)
            ring_plot_magnitudes.append(plot_magnitudes_set)
        plot_magnitudes_by_rsv.append(ring_plot_magnitudes)

    # Calculate plot point coordinates
    line_points_by_rvs = [ ] # RVS = ring, version, segment
    discrete_segs_by_rvs = [ ]
    for ring_id in range(num_rings):
        ring_base_radius = (ring_id + 2) * division_size
        # Continuous lines
        if settings['ring_types'][ring_id] == 'C':
            line_points = [ [ ] for _ in range(num_versions) ]
            # Add zero-point
            for version_id in range(num_versions):
                point = _coords_from_angle(center, angle_delta*0.5, ring_base_radius, gap=gap)
                line_points[version_id].append(point)
            for segment_id in range(num_segments):
                base_angle = angle_delta * (segment_id+0.5)
                for version_id, dp in enumerate(dp_set):
                    plot_magnitude = plot_magnitudes_by_rsv[ring_id][segment_id][version_id]
                    if plot_magnitude is None:
                        point = _coords_from_angle(center, base_angle, ring_base_radius, gap=gap)
                    else:
                        plot_radius = ring_base_radius + division_size * plot_magnitude
                        point = _coords_from_angle(center, base_angle, plot_radius, gap=gap)
                    line_points[version_id].append(point)
            line_points_by_rvs.append(line_points)
            discrete_segs_by_rvs.append(None)
        # Discrete segments
        elif settings['ring_types'][ring_id] == 'D':
            ring_discrete_segs = [ [ ] for _ in range(num_versions) ]
            for segment_id, dp_set in enumerate(datapoints):
                base_angle = angle_delta * (segment_id+0.5)
                for version_id, dp in enumerate(dp_set):
                    segment_length = sizes[0]
                    segment_color = settings['discrete_color_sequence'][-1]
                    segment_opacity = 1
                    if dp is not None and dp['discrete'] is not None and dp['discrete'][ring_id] is not None:
                        category_id = dp['discrete'][ring_id]
                        segment_color = settings['discrete_color_sequence'][category_id]
                    if segment_color == settings['discrete_color_sequence'][-1]:
                        segment_opacity = 0.5
                    segment_points = (_coords_from_angle(center, angle_delta*(segment_id), ring_base_radius - segment_length, gap=gap),
                                      _coords_from_angle(center, angle_delta*(segment_id), ring_base_radius + segment_length, gap=gap),
                                      _coords_from_angle(center, angle_delta*(segment_id+1), ring_base_radius + segment_length, gap=gap),
                                      _coords_from_angle(center, angle_delta*(segment_id+1), ring_base_radius - segment_length, gap=gap))
                    ring_discrete_segs[version_id].append((segment_points, segment_color, segment_opacity))
            line_points_by_rvs.append(None)
            discrete_segs_by_rvs.append(ring_discrete_segs)

    # Draw plot points
    for ring_id in range(num_rings):
        baseline_circle_points = [ ]
        for i in range(settings['baseline_point_resolution'] + 1):
            point_angle = (settings['baseline_point_resolution'] - i) * (2*pi - gap) / float(settings['baseline_point_resolution'])
            point_radius = (ring_id + 2) * division_size
            baseline_circle_points.append(_coords_from_angle(center, point_angle, point_radius, gap=gap))
        # Continuous lines
        if settings['ring_types'][ring_id] == 'C':
            line_points = line_points_by_rvs[ring_id][-1] + baseline_circle_points
            ring_line = dwg.polyline(line_points,
                                     stroke=settings['ring_color_sequence'][ring_id],
                                     stroke_width=2,
                                     stroke_opacity=1,
                                     fill=settings['ring_color_sequence'][ring_id],
                                     fill_opacity=0.2,
                                     id=settings['svg_id'] + '-contline-' + str(ring_id))
            for version_id in range(num_versions):
                line_points = line_points_by_rvs[ring_id][version_id] + baseline_circle_points
                points_string = ' '.join([ ','.join([ str(c) for c in point ]) for point in line_points ])
                animation = Animate(values=None,
                                    dur=str(settings['animation_length']) + 'ms',
                                    begin='indefinite',
                                    fill='freeze',
                                    attributeName='points',
                                    to=points_string,
                                    id=settings['svg_id'] + '-animation-' + str(ring_id) + '-' + str(version_id))
                ring_line.add(animation)
            dwg.add(ring_line)
        # Discrete segments
        elif settings['ring_types'][ring_id] == 'D':
            for version_id in range(num_versions):
                group_opacity = 1 if version_id == num_versions-1 else 0
                segment_group = dwg.g(id=settings['svg_id'] + '-discrete-' + str(ring_id) + '-' + str(version_id), opacity=group_opacity)
                segment_point_groups = discrete_segs_by_rvs[ring_id][version_id]
                for segment_points, segment_color, segment_opacity in segment_point_groups:
                    segment_group.add(dwg.polyline(segment_points,
                                                   stroke_width=0,
                                                   stroke_opacity=0,
                                                   fill=segment_color,
                                                   fill_opacity=segment_opacity))
                dwg.add(segment_group)

    # Draw missing-data shade
    if settings['missing_data_shade_enabled']:
        for version_id in range(num_versions):
            group_opacity = 1 if version_id == num_versions-1 else 0
            shade_group = dwg.g(id=settings['svg_id'] + '-shade-' + str(version_id), opacity=group_opacity)
            for dp_id, dp_set in enumerate(datapoints):
                if dp_set[version_id] is None:
                    shade_group.add(dwg.polygon([ center,
                                                  _coords_from_angle(center, angle_delta*dp_id, full_radius+5, gap=gap),
                                                  _coords_from_angle(center, angle_delta*(dp_id+1), full_radius+5, gap=gap) ],
                                                  stroke_opacity=0,
                                                  fill=settings['missing_data_shade_color'],
                                                  fill_opacity=1))
                    """
                    shade_group.add(dwg.circle(r=sizes[0]//2,
                                               center=_coords_from_angle(center, angle_delta*(dp_id+0.5), full_radius-1.5*sizes[1], gap=gap),
                                               fill=COLORS['BLACK'],
                                               fill_opacity=1,
                                               stroke_opacity=0))
                    """
            dwg.add(shade_group)

    # Draw outer markers
    if settings['apply_markers']:
        any_markers = False
        for version_id in range(num_versions):
            group_opacity = 1 if version_id == num_versions-1 else 0
            markers_group = dwg.g(id=settings['svg_id'] + '-markers-' + str(version_id), opacity=group_opacity)
            for dp_id, dp_set in enumerate(datapoints):
                if dp_set[version_id] is None:
                    continue
                if 'marker' in dp_set[version_id] and dp_set[version_id]['marker'] == True:
                    any_markers = True
                    markers_group.add(dwg.line(_coords_from_angle(center, angle_delta * (dp_id+0.2), full_radius-sizes[0]*0.1, gap=gap),
                                               _coords_from_angle(center, angle_delta * (dp_id+0.8), full_radius-sizes[0]*0.9, gap=gap),
                                               stroke=settings['marker_color'],
                                               stroke_width=3,
                                               stroke_opacity=1))
                    markers_group.add(dwg.line(_coords_from_angle(center, angle_delta * (dp_id+0.2), full_radius-sizes[0]*0.9, gap=gap),
                                               _coords_from_angle(center, angle_delta * (dp_id+0.8), full_radius-sizes[0]*0.1, gap=gap),
                                               stroke=settings['marker_color'],
                                               stroke_width=3,
                                               stroke_opacity=1))
            dwg.add(markers_group)
        if any_markers:
            dwg.add(dwg.text(text=settings['marker_label'],
                             insert=_coords_from_angle(center, 0, full_radius-sizes[0]*0.1),
                             font_size=sizes[1],
                             font_family='Arial',
                             text_anchor='middle',
                             alignment_baseline='central'))

    # Draw outer rings
    dwg.add(dwg.circle(r=full_radius-2*sizes[1],
                       center=center,
                       fill_opacity=0,
                       stroke=COLORS['BLACK'],
                       stroke_width=1,
                       stroke_opacity=0.5))
    dwg.add(dwg.circle(r=full_radius-1*sizes[1],
                       center=center,
                       fill_opacity=0,
                       stroke=COLORS['BLACK'],
                       stroke_width=1,
                       stroke_opacity=0.5))
    for i in range(num_segments+1):
        dwg.add(dwg.line(_coords_from_angle(center, angle_delta*i, full_radius - 2*sizes[1], gap=gap),
                         _coords_from_angle(center, angle_delta*i, full_radius - 1*sizes[1], gap=gap),
                         stroke=COLORS['BLACK'],
                         stroke_width=1,
                         stroke_opacity=0.5))

    # Draw segment selector
    if settings['segment_selector_enabled']:
        center_point = angle_delta*0.5
        selector_points = (_coords_from_angle(center, center_point, full_radius-1.5*sizes[1], gap=gap),
                           _coords_from_angle(center, center_point-0.02, full_radius-0.5*sizes[1], gap=gap),
                           _coords_from_angle(center, center_point-0.02, full_radius+5, gap=gap),
                           _coords_from_angle(center, center_point+0.02, full_radius+5, gap=gap),
                           _coords_from_angle(center, center_point+0.02, full_radius-0.5*sizes[1], gap=gap))
        dwg.add(dwg.polygon(selector_points,
                            stroke=COLORS['BLACK'],
                            stroke_width=2,
                            stroke_opacity=1,
                            fill=COLORS['D_GREY'],
                            fill_opacity=0.2,
                            id=settings['svg_id'] + '-selector'))

    # Draw interaction segments
    if settings['interaction_segments_enabled']:
        for i in range(num_segments):
            dwg.add(dwg.polygon([ center,
                                  _coords_from_angle(center, angle_delta*i, full_radius+5, gap=gap),
                                  _coords_from_angle(center, angle_delta*(i+1), full_radius+5, gap=gap) ],
                                  stroke=COLORS['BLACK'],
                                  stroke_width=1,
                                  stroke_opacity=0,
                                  fill=COLORS['L_GREY'],
                                  fill_opacity=0,
                                  onmousedown=settings['segment_function_name'] + '(1, ' + str(i) +');',
                                  onmouseover=settings['segment_function_name'] + '(2, ' + str(i) +');',
                                  onmouseup=settings['segment_function_name'] + '(3, ' + str(i) +');',
                                  id=settings['svg_id'] + '-intseg-' + str(i)))
        dwg.add(dwg.circle(r=1.5*division_size,
                           center=center,
                           fill=COLORS['WHITE'],
                           fill_opacity=1,
                           stroke_opacity=0))

    # Draw center text
    if settings['center_text_1'] is not None:
        dwg.add(dwg.text(text=settings['center_text_1'],
                         insert=(center[0], center[1]-1.5*sizes[1]),
                         font_size=1.5*sizes[1],
                         font_family='Arial',
                         font_weight='bold',
                         text_anchor='middle',
                         alignment_baseline='central'))
    if settings['center_text_2'] is not None:
        dwg.add(dwg.text(text=settings['center_text_2'],
                         insert=(center[0], center[1]+1*sizes[1]),
                         font_size=sizes[1],
                         font_family='Arial',
                         text_anchor='middle',
                         alignment_baseline='central'))
    if settings['center_text_3'] is not None:
        dwg.add(dwg.text(text=settings['center_text_3'],
                         insert=(center[0], center[1]+3*sizes[1]),
                         font_size=sizes[1],
                         font_family='Arial',
                         text_anchor='middle',
                         alignment_baseline='central',
                         fill=COLORS['L_GREY']))
    return dwg.tostring()


def radar(axis_names=('A', 'B', 'C', 'D', 'E'), settings={ }):
    # Check settings
    for setting_name, setting_default in DEFAULT_SETTINGS_RADAR.items():
        if setting_name not in settings.keys():
            settings[setting_name] = setting_default

    # Useful calculations
    num_metrics = len(axis_names)
    center = (settings['canvas_size'][0]//2, settings['canvas_size'][1]//2)
    sizes = (min(settings['canvas_size'])//50, min(settings['canvas_size'])//40)
    angle_delta = 2*pi / float(num_metrics)
    axis_radius = min(settings['canvas_size'])//2 - min(settings['canvas_size'])//20

    # Initialise drawing
    dwg = svgwrite.Drawing(profile='full')

    # Set viewbox attribute for scaling
    dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in settings['canvas_size'] ])
    dwg.attribs['width'] = '100%'
    dwg.attribs['height'] = '100%'

    # Set HTML attributes for interaction
    # If the user doesn't set an ID for the chart, then it shouldn't have an ID defined in the XML. But, the svg_id variable still needs to be
    # defined for later, to make sure there's something to prepend to the IDs of any child elements that *do* need IDs.
    if settings['svg_id'] is None:
        settings['svg_id'] = 'radar'
    else:
        dwg.attribs['id'] = settings['svg_id']
    if settings['svg_hidden']:
        dwg.attribs['style'] = 'display: none;'

    # Axes
    for i in range(num_metrics):
        dwg.add(dwg.line(center,
                         _coords_from_angle(center, angle_delta*i, axis_radius),
                         stroke=COLORS['L_GREY'],
                         stroke_width=2))
        for j in range(10):
            dwg.add(dwg.line(_coords_from_angle(center, angle_delta*i, (j+1)*axis_radius//10),
                             _coords_from_angle(center, angle_delta*(i+1), (j+1)*axis_radius//10),
                             stroke=COLORS['L_GREY'],
                             stroke_width=1))

    # Axis labels
    for i in range(num_metrics):
        l_width, l_height = sizes[1]/float(2) * len(axis_names[i]), sizes[1]
        vert_adj = l_height if pi/2 < angle_delta*i < 3*pi/2 else -l_height if 3*pi/2 < angle_delta*i or angle_delta*i >= 0 else 0
        horiz_adj = l_width/2 if 0 < angle_delta*i < pi else -l_width/2 if pi < angle_delta*i < 2*pi else 0
        dwg.add(dwg.text(text=axis_names[i],
                         insert=_coords_from_angle(center, angle_delta*i, axis_radius, adj=(horiz_adj, vert_adj)),
                         font_size=sizes[1],
                         font_family='Arial',
                         text_anchor='middle',
                         alignment_baseline='central'))
    for i in range(10):
        label = str(10*(i+1))
        dwg.add(dwg.rect(insert=(center[0] - sizes[0], center[1] - (i+1)*axis_radius//10 - sizes[0]//2),
                         size=(sizes[0]*2, sizes[0]*1.5),
                         fill=COLORS['WHITE']))
        dwg.add(dwg.text(text=label,
                         insert=(center[0] - sizes[0]//1.5 - sizes[0]//3*(len(label)-2), center[1] - (i+1)*axis_radius//10 + sizes[0]//2),
                         font_size=sizes[0],
                         font_family='Arial',
                         fill=COLORS['D_GREY']))

    # Dots, labels, and connecting polygon (initialised hidden)
    dwg.add(dwg.polygon(points=[ _coords_from_angle(center, angle_delta*i, axis_radius/2) for i in range(len(axis_names)) ],
                        fill=COLORS['WHITE'],
                        fill_opacity=0,
                        stroke=COLORS['BLACK'],
                        stroke_width=1,
                        stroke_opacity=0,
                        id='connecting-polygon'))
    for i in range(len(axis_names)):
        dwg.add(dwg.circle(r=sizes[0],
                           center=center,
                           fill=COLORS['L_GREY'],
                           fill_opacity=0,
                           stroke=COLORS['BLACK'],
                           stroke_width=1,
                           stroke_opacity=0,
                           id='metric-point-' + str(i),
                           onclick='handleRadarChart(1, ' + str(i) +');',
                           onmouseover='handleRadarChart(2, ' + str(i) +');',
                           onmouseout='handleRadarChart(3, ' + str(i) +');'))
    floating_label = dwg.g(id='floating-label-group',
                           opacity=0)
    floating_label.add(dwg.circle(r=sizes[1]*1.66,
                                  center=(0, 0),
                                  stroke=COLORS['D_GREY'],
                                  stroke_width=1,
                                  fill=COLORS['WHITE'],
                                  fill_opacity=0.9))
    floating_label.add(dwg.text('',
                                insert=(0, -sizes[0]*0.66),
                                font_size=sizes[0],
                                font_family='Arial',
                                font_weight='bold',
                                fill=COLORS['BLACK'],
                                fill_opacity=1,
                                text_anchor='middle',
                                alignment_baseline='central',
                                id='floating-label-text-top'))
    floating_label.add(dwg.text('',
                                insert=(0, +sizes[0]*0.66),
                                font_size=sizes[0],
                                font_family='Arial',
                                fill=COLORS['BLACK'],
                                fill_opacity=1,
                                text_anchor='middle',
                                alignment_baseline='central',
                                id='floating-label-text-bottom'))
    dwg.add(floating_label)

    return dwg.tostring()


def grid(settings={ }):
    # Check arguments
    for setting_name, setting_default in DEFAULT_SETTINGS_GRID.items():
        if setting_name not in settings.keys():
            settings[setting_name] = setting_default

    # Useful calculations
    # Bounds defined as (x0, y0, x1, y1)
    sizes = (min(settings['canvas_size'])//50, min(settings['canvas_size'])//40)
    margin = 2*sizes[0]
    bar_width = 15*sizes[0]
    tickbox_1_bounds = (margin, margin, settings['canvas_size'][0]*0.45, 15*sizes[0])
    tickbox_2_bounds = (settings['canvas_size'][0]*0.55, margin, settings['canvas_size'][0]-margin, 15*sizes[0])
    divider_line_y = 22.5*sizes[0]
    bar_charts_bounds = (margin+3*sizes[0], 28*sizes[0], settings['canvas_size'][0]-margin, settings['canvas_size'][1]-10*sizes[0])
    bar_1_x = bar_charts_bounds[0] + (bar_charts_bounds[2]-bar_charts_bounds[0])*1/4
    bar_2_x = bar_charts_bounds[0] + (bar_charts_bounds[2]-bar_charts_bounds[0])*3/4

    # Initialise drawing
    dwg = svgwrite.Drawing(profile='full')

    # Set viewbox attribute for scaling
    dwg.attribs['viewBox'] = '0 0 ' + ' '.join([ str(x) for x in settings['canvas_size'] ])
    dwg.attribs['width'] = '95%'
    dwg.attribs['height'] = '95%'

    # Set HTML attributes for interaction
    # If the user doesn't set an ID for the chart, then it shouldn't have an ID defined in the XML. But, the svg_id variable still needs to be
    # defined for later, to make sure there's something to prepend to the IDs of any child elements that *do* need IDs.
    if settings['svg_id'] is None:
        settings['svg_id'] = 'grid'
    else:
        dwg.attribs['id'] = settings['svg_id']
    if settings['svg_hidden']:
        dwg.attribs['style'] = 'display: none;'

    # Checkboxes
    dwg.add(dwg.polygon(points=[ (tickbox_1_bounds[0], tickbox_1_bounds[1]),
                                 (tickbox_1_bounds[2], tickbox_1_bounds[1]),
                                 (tickbox_1_bounds[2], tickbox_1_bounds[3]),
                                 (tickbox_1_bounds[0], tickbox_1_bounds[3]) ],
                        fill=COLORS['VL_GREY'],
                        fill_opacity=0.8,
                        stroke=COLORS['BLACK'],
                        stroke_width=2,
                        stroke_opacity=1,
                        id='checkbox-1'))
    dwg.add(dwg.text('',
                      insert=((tickbox_1_bounds[0]+tickbox_1_bounds[2])/2, (tickbox_1_bounds[1]+tickbox_1_bounds[3])/2),
                      font_size=2*sizes[1],
                      font_family='Arial',
                      font_weight='bold',
                      fill=COLORS['BLACK'],
                      fill_opacity=1,
                      text_anchor='middle',
                      alignment_baseline='central',
                      id='checkbox-1-text'))
    dwg.add(dwg.text(settings['box_1_label'],
                      insert=((tickbox_1_bounds[0]+tickbox_1_bounds[2])/2, tickbox_1_bounds[3]+25),
                      font_size=1.8*sizes[1],
                      font_family='Arial',
                      fill=COLORS['BLACK'],
                      fill_opacity=1,
                      text_anchor='middle',
                      alignment_baseline='central'))
    dwg.add(dwg.polygon(points=[ (tickbox_2_bounds[0], tickbox_2_bounds[1]),
                                 (tickbox_2_bounds[2], tickbox_2_bounds[1]),
                                 (tickbox_2_bounds[2], tickbox_2_bounds[3]),
                                 (tickbox_2_bounds[0], tickbox_2_bounds[3]) ],
                        fill=COLORS['VL_GREY'],
                        fill_opacity=0.8,
                        stroke=COLORS['BLACK'],
                        stroke_width=2,
                        stroke_opacity=1,
                        id='checkbox-2'))
    dwg.add(dwg.text('',
                      insert=((tickbox_2_bounds[0]+tickbox_2_bounds[2])/2, (tickbox_2_bounds[1]+tickbox_2_bounds[3])/2),
                      font_size=2*sizes[1],
                      font_family='Arial',
                      font_weight='bold',
                      fill=COLORS['BLACK'],
                      fill_opacity=1,
                      text_anchor='middle',
                      alignment_baseline='central',
                      id='checkbox-2-text'))
    dwg.add(dwg.text(settings['box_2_label'],
                      insert=((tickbox_2_bounds[0]+tickbox_2_bounds[2])/2, tickbox_2_bounds[3]+25),
                      font_size=1.8*sizes[1],
                      font_family='Arial',
                      fill=COLORS['BLACK'],
                      fill_opacity=1,
                      text_anchor='middle',
                      alignment_baseline='central'))

    # Divider line
    dwg.add(dwg.line((margin, divider_line_y),
                     (settings['canvas_size'][0]-margin, divider_line_y),
                     stroke=COLORS['L_GREY'],
                     stroke_width=2))

    # Bars
    dwg.add(dwg.polygon(points=[ (bar_charts_bounds[0], bar_charts_bounds[1]),
                                 (bar_charts_bounds[2], bar_charts_bounds[1]),
                                 (bar_charts_bounds[2], bar_charts_bounds[3]),
                                 (bar_charts_bounds[0], bar_charts_bounds[3]) ],
                        fill=COLORS['WHITE'],
                        fill_opacity=0,
                        stroke=COLORS['BLACK'],
                        stroke_width=2,
                        stroke_opacity=1,
                        id='bar-charts-container'))
    dwg.add(dwg.text(settings['bar_label'],
                      insert=(settings['canvas_size'][0]/2, bar_charts_bounds[3]+6*sizes[0]),
                      font_size=1.8*sizes[1],
                      font_family='Arial',
                      fill=COLORS['BLACK'],
                      fill_opacity=1,
                      text_anchor='middle',
                      alignment_baseline='central'))
    dwg.add(dwg.text(settings['bar_1_label'],
                      insert=(bar_1_x, bar_charts_bounds[3]+20),
                      font_size=1.8*sizes[1],
                      font_family='Arial',
                      fill=COLORS['BLACK'],
                      fill_opacity=1,
                      text_anchor='middle',
                      alignment_baseline='central'))
    dwg.add(dwg.text(settings['bar_2_label'],
                      insert=(bar_2_x, bar_charts_bounds[3]+20),
                      font_size=1.8*sizes[1],
                      font_family='Arial',
                      fill=COLORS['BLACK'],
                      fill_opacity=1,
                      text_anchor='middle',
                      alignment_baseline='central'))
    # Bar axis
    for i in range(10+1):
        height = bar_charts_bounds[1]+i*(bar_charts_bounds[3]-bar_charts_bounds[1])/10
        dwg.add(dwg.line((bar_charts_bounds[0]-5, height), (bar_charts_bounds[0]+5, height),
                         stroke=COLORS['BLACK'],
                         stroke_width=2,
                         stroke_opacity=1))
        dwg.add(dwg.text(str(100-i*10),
                         insert=(bar_charts_bounds[0]-sizes[0], height),
                         font_size=1.5*sizes[1],
                         font_family='Arial',
                         fill=COLORS['BLACK'],
                         fill_opacity=1,
                         text_anchor='end',
                         alignment_baseline='central'))
    # Bar 1
    dwg.add(dwg.polygon(points=[ (bar_1_x-bar_width//2, bar_charts_bounds[3]),
                                 (bar_1_x-bar_width//2, bar_charts_bounds[1]),
                                 (bar_1_x+bar_width//2, bar_charts_bounds[1]),
                                 (bar_1_x+bar_width//2, bar_charts_bounds[3]) ],
                        fill=COLORS['VL_GREY'],
                        fill_opacity=0.5,
                        stroke=COLORS['BLACK'],
                        stroke_width=2,
                        stroke_opacity=1))
    # Box plot 1
    box_plot_group_1 = dwg.g(id='boxplot-1', opacity=0)
    box_plot_group_1.add(dwg.polygon(points=[ (bar_1_x-bar_width//2, bar_charts_bounds[3]),
                                              (bar_1_x-bar_width//2, bar_charts_bounds[1]),
                                              (bar_1_x+bar_width//2, bar_charts_bounds[1]),
                                              (bar_1_x+bar_width//2, bar_charts_bounds[3]) ],
                                    fill=COLORS['WHITE'],
                                    fill_opacity=1,
                                    stroke=COLORS['BLACK'],
                                    stroke_width=2,
                                    stroke_opacity=1))
    box_plot_group_1.add(dwg.polygon(points=[ (bar_1_x-bar_width//2, bar_charts_bounds[1]+80),
                                              (bar_1_x-bar_width//2, bar_charts_bounds[3]-80),
                                              (bar_1_x+bar_width//2, bar_charts_bounds[3]-80),
                                              (bar_1_x+bar_width//2, bar_charts_bounds[1]+80) ],
                                    fill='url(#gradient-1)',
                                    fill_opacity=0.8,
                                    stroke=COLORS['BLACK'],
                                    stroke_width=2,
                                    stroke_opacity=0.5,
                                    id='boxplot-1-box'))
    gradient = LinearGradient(start=(0, 0), end=(0,1), id='gradient-1')
    gradient.add_stop_color(offset='0%', color=COLORS['BAR_GREEN'])
    gradient.add_stop_color(offset='50%', color=COLORS['BAR_ORANGE'])
    gradient.add_stop_color(offset='100%', color=COLORS['BAR_RED'])
    dwg.defs.add(gradient)
    box_plot_group_1.add(dwg.line((bar_1_x-bar_width//2, bar_charts_bounds[1]+200),
                                  (bar_1_x+bar_width//2, bar_charts_bounds[1]+200),
                                  stroke=COLORS['BLACK'],
                                  stroke_width=2,
                                  stroke_opacity=0.5,
                                  stroke_dasharray=2,
                                  id='boxplot-1-line-high'))
    box_plot_group_1.add(dwg.line((bar_1_x-bar_width//2, (bar_charts_bounds[1]+bar_charts_bounds[3])//2),
                                  (bar_1_x+bar_width//2, (bar_charts_bounds[1]+bar_charts_bounds[3])//2),
                                  stroke=COLORS['BLACK'],
                                  stroke_width=3,
                                  stroke_opacity=0.8,
                                  stroke_dasharray=5,
                                  id='boxplot-1-line-mid'))
    box_plot_group_1.add(dwg.line((bar_1_x-bar_width//2, bar_charts_bounds[3]-200),
                                  (bar_1_x+bar_width//2, bar_charts_bounds[3]-200),
                                  stroke=COLORS['BLACK'],
                                  stroke_width=2,
                                  stroke_opacity=0.5,
                                  stroke_dasharray=2,
                                  id='boxplot-1-line-low'))
    box_plot_group_1.add(dwg.line((bar_1_x-bar_width//2, bar_charts_bounds[3]),
                                  (bar_1_x+bar_width//2, bar_charts_bounds[3]),
                                  fill_opacity=0,
                                  stroke=COLORS['BLACK'],
                                  stroke_width=4,
                                  stroke_opacity=1,
                                  id='bar-1-mainline'))
    box_plot_group_1.add(dwg.text('',
                                  insert=(bar_1_x, bar_charts_bounds[3]),
                                  font_size=2*sizes[1],
                                  font_family='Arial',
                                  font_weight='bold',
                                  fill=COLORS['BLACK'],
                                  fill_opacity=1,
                                  text_anchor='middle',
                                  alignment_baseline='central',
                                  id='bar-1-label'))
    dwg.add(box_plot_group_1)
    # Bar 2
    dwg.add(dwg.polygon(points=[ (bar_2_x-bar_width//2, bar_charts_bounds[3]),
                                 (bar_2_x-bar_width//2, bar_charts_bounds[1]),
                                 (bar_2_x+bar_width//2, bar_charts_bounds[1]),
                                 (bar_2_x+bar_width//2, bar_charts_bounds[3]) ],
                        fill=COLORS['VL_GREY'],
                        fill_opacity=0.5,
                        stroke=COLORS['BLACK'],
                        stroke_width=2,
                        stroke_opacity=1))
    # Box plot 2
    box_plot_group_2 = dwg.g(id='boxplot-2', opacity=0)
    box_plot_group_2.add(dwg.polygon(points=[ (bar_2_x-bar_width//2, bar_charts_bounds[3]),
                                              (bar_2_x-bar_width//2, bar_charts_bounds[1]),
                                              (bar_2_x+bar_width//2, bar_charts_bounds[1]),
                                              (bar_2_x+bar_width//2, bar_charts_bounds[3]) ],
                                    fill=COLORS['WHITE'],
                                    fill_opacity=1,
                                    stroke=COLORS['BLACK'],
                                    stroke_width=2,
                                    stroke_opacity=1))
    box_plot_group_2.add(dwg.polygon(points=[ (bar_2_x-bar_width//2, bar_charts_bounds[1]+80),
                                              (bar_2_x-bar_width//2, bar_charts_bounds[3]-80),
                                              (bar_2_x+bar_width//2, bar_charts_bounds[3]-80),
                                              (bar_2_x+bar_width//2, bar_charts_bounds[1]+80) ],
                                    fill='url(#gradient-2)',
                                    fill_opacity=0.8,
                                    stroke=COLORS['BLACK'],
                                    stroke_width=2,
                                    stroke_opacity=0.5,
                                    id='boxplot-2-box'))
    gradient = LinearGradient(start=(0, 0), end=(0,1), id='gradient-2')
    gradient.add_stop_color(offset='0%', color=COLORS['BAR_GREEN'])
    gradient.add_stop_color(offset='50%', color=COLORS['BAR_ORANGE'])
    gradient.add_stop_color(offset='100%', color=COLORS['BAR_RED'])
    dwg.defs.add(gradient)
    box_plot_group_2.add(dwg.line((bar_2_x-bar_width//2, bar_charts_bounds[1]+200),
                                  (bar_2_x+bar_width//2, bar_charts_bounds[1]+200),
                                  stroke=COLORS['BLACK'],
                                  stroke_width=2,
                                  stroke_opacity=0.5,
                                  stroke_dasharray=2,
                                  id='boxplot-2-line-high'))
    box_plot_group_2.add(dwg.line((bar_2_x-bar_width//2, (bar_charts_bounds[1]+bar_charts_bounds[3])//2),
                                  (bar_2_x+bar_width//2, (bar_charts_bounds[1]+bar_charts_bounds[3])//2),
                                  stroke=COLORS['BLACK'],
                                  stroke_width=3,
                                  stroke_opacity=0.8,
                                  stroke_dasharray=5,
                                  id='boxplot-2-line-mid'))
    box_plot_group_2.add(dwg.line((bar_2_x-bar_width//2, bar_charts_bounds[3]-200),
                                  (bar_2_x+bar_width//2, bar_charts_bounds[3]-200),
                                  stroke=COLORS['BLACK'],
                                  stroke_width=2,
                                  stroke_opacity=0.5,
                                  stroke_dasharray=2,
                                  id='boxplot-2-line-low'))
    box_plot_group_2.add(dwg.line((bar_2_x-bar_width//2, bar_charts_bounds[3]),
                                  (bar_2_x+bar_width//2, bar_charts_bounds[3]),
                                  fill_opacity=0,
                                  stroke=COLORS['BLACK'],
                                  stroke_width=4,
                                  stroke_opacity=1,
                                  id='bar-2-mainline'))
    box_plot_group_2.add(dwg.text('',
                                  insert=(bar_2_x, bar_charts_bounds[3]),
                                  font_size=2*sizes[1],
                                  font_family='Arial',
                                  font_weight='bold',
                                  fill=COLORS['BLACK'],
                                  fill_opacity=1,
                                  text_anchor='middle',
                                  alignment_baseline='central',
                                  id='bar-2-label'))
    dwg.add(box_plot_group_2)

    return dwg.tostring()
