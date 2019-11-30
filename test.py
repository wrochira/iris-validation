#!/usr/bin/env python3

from multimetric_report import model
from multimetric_report import report

protein = model.generate(chains=3, chain_length=500, chain_flaws=5)
report.generate(protein, sidepanel_height=50, output_dir='../REPORT')
