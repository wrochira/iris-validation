"""
Copyright 2020 William Rochira at the University of York
"""

from clipper_tools.metrics import generate_metrics_model

from .stats import chart_data_from_models
from .report import report_from_data


def generate_report(model_paths, reflections_paths, output_path, mode='', js_old=False):
    if len(model_paths) != len(reflections_paths):
        print('ERROR: path arguments should both be of the same length')
        return
    if len(model_paths) not in (1, 2):
        print('ERROR: path arguments should be iterables of length 1 or 2')
        return
    metrics_models = [ ]
    for model_path, reflections_path in zip(model_paths, reflections_paths):
        metrics_model = generate_metrics_model(model_path, reflections_path)
        metrics_models.append(metrics_model)
    chart_data = chart_data_from_models(*reversed(metrics_models))
    if chart_data is None:
        print('ERROR: failed to generate chart data from MetricsModel')
        return
    report_from_data(chart_data, output_path, mode, js_old)
