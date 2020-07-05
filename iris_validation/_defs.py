"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

NO_CONTEXT_MESSAGE = 'NO CONTEXT'
SC_INCOMPLETE_STRING = 'INCOMPLETE SIDECHAIN'

RESOLUTION_BIN_NAMES = ('<10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '>90', 'All')

METRIC_NAMES = ('Ramachandran Score', 'Rotamer Score', 'Avg B-factor', 'Max B-factor', 'Std B-factor', 'Residue Fit', 'Mainchain Fit', 'Sidechain Fit')
METRIC_POLARITIES = (+1, -1, -1, -1, -1, -1, -1, -1)
METRIC_SHORTNAMES = ('Rama', 'Rota', 'Avg B', 'Max B', 'Std B', 'Res. Fit', 'MC Fit', 'SC Fit')
METRIC_DISPLAY_TYPES = ('D', 'D', 'C', 'C', 'C', 'C', 'C', 'C')

REPORT_METRIC_IDS = (0, 1, 2, 3, 6, 7)
REPORT_RESIDUE_VIEW = 'Grid'
