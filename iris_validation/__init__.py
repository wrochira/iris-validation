"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

from iris_validation._defs import *
from iris_validation import metrics, interface


def generate_report(latest_model_path,
                    previous_model_path=None,
                    latest_reflections_path=None,
                    previous_reflections_path=None,
                    output_dir='./iris_report_output/',
                    mode=''):
    
    latest_metrics_model = metrics.generate_metrics_model(latest_model_path, latest_reflections_path)
    previous_metrics_model = None
    if previous_model_path is not None:
        previous_metrics_model = metrics.generate_metrics_model(previous_model_path, previous_reflections_path)

    interface.build_report(latest_metrics_model, previous_metrics_model, output_dir, mode)
