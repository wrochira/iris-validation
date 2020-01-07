from random import random, shuffle

import numpy as np
from scipy.stats import norm


RESIDUES = [ 'A', 'R', 'N', 'D', 'B', 'C', 'E', 'Q', 'Z', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' ]


class ReportResidue(object):
	def __init__(self, code, metrics):
		super(ReportResidue, self).__init__()
		self.code = code
		self.metrics = metrics


class ReportChain(object):
	def __init__(self, chain_id, residues):
		super(ReportChain, self).__init__()
		self.chain_id = chain_id
		self.residues = residues
		self.length = len(residues)


class ReportProtein(object):
	def __init__(self, chains):
		super(ReportProtein, self).__init__()
		self.chains = chains
		self.all_residues = [ residue for residues in [ [ residue for residue in chain.residues ] for chain in chains ] for residue in residues ]
		try:
			self.metric_names = [ metric['name'] for metric in self.all_residues[0].metrics ]
			self.metric_shortnames = [ metric['shortname'] for metric in self.all_residues[0].metrics ]
		except IndexError:
			raise Exception('Protein is formatted incorrectly')
			exit(1)

		self.metric_stats_abs = { }
		self.metric_stats_rep = { }
		self.metric_quartiles_abs = { }
		self.metric_quartiles_rep = { }
		self.metric_sqdiff_minmax_abs = { }
		self.metric_sqdiff_minmax_rep = { }

		for i in range(len(self.metric_names)):
			all_values_abs = [ residue.metrics[i]['absolute_value'] for residue in self.all_residues if residue.metrics[i]['absolute_value'] is not None ]
			all_values_rep = [ residue.metrics[i]['represented_value'] for residue in self.all_residues if residue.metrics[i]['represented_value'] is not None ]
			if len(all_values_abs) == 0 or len(all_values_rep) == 0:
				raise Exception('One or more metrics have no values specified')
			mean_abs, sdev_abs = sum(all_values_abs)/len(all_values_abs), np.std(all_values_abs)
			mean_rep, sdev_rep = sum(all_values_rep)/len(all_values_rep), np.std(all_values_rep)
			self.metric_stats_abs[i] = (mean_abs, sdev_abs)
			self.metric_stats_rep[i] = (mean_rep, sdev_rep)
			self.metric_quartiles_abs[i] = [ np.percentile(all_values_abs, q) for q in (0, 25, 50, 75, 100) ]
			self.metric_quartiles_rep[i] = [ np.percentile(all_values_rep, q) for q in (0, 25, 50, 75, 100) ]
			metric_sqdiff_polarities_abs = (1 if self.metric_quartiles_abs[i][0] > mean_abs else -1, 1 if self.metric_quartiles_abs[i][-1] > mean_abs else -1)
			metric_sqdiff_polarities_rep = (1 if self.metric_quartiles_rep[i][0] > mean_rep else -1, 1 if self.metric_quartiles_rep[i][-1] > mean_rep else -1)
			self.metric_sqdiff_minmax_abs[i] = (metric_sqdiff_polarities_abs[0]*((self.metric_quartiles_abs[i][0]-mean_abs)/sdev_abs)**2, metric_sqdiff_polarities_abs[1]*((self.metric_quartiles_abs[i][-1]-mean_abs)/sdev_abs)**2)
			self.metric_sqdiff_minmax_rep[i] = (metric_sqdiff_polarities_rep[0]*((self.metric_quartiles_rep[i][0]-mean_rep)/sdev_rep)**2, metric_sqdiff_polarities_rep[1]*((self.metric_quartiles_rep[i][-1]-mean_rep)/sdev_rep)**2)


class ReportProteinGroup(object):
	def __init__(self, proteins):
		super(ReportProteinGroup, self).__init__()
		self.proteins = proteins
		self.all_residues = [ residue for all_residues in [ protein.all_residues for protein in proteins ] for residue in all_residues ]
		try:
			self.metric_names = [ metric['name'] for metric in self.all_residues[0].metrics ]
			self.metric_shortnames = [ metric['shortname'] for metric in self.all_residues[0].metrics ]
		except IndexError:
			raise Exception('Protein is formatted incorrectly')
		if [ len(protein.chains) for protein in proteins ].count(len(protein.chains)) != len(proteins):
			raise Exception('Model ghosting may be incorrect for this model series due to differening chain counts')

		self.metric_stats_abs = { }
		self.metric_stats_rep = { }
		self.metric_quartiles_abs = { }
		self.metric_quartiles_rep = { }
		self.metric_sqdiff_minmax_abs = { }
		self.metric_sqdiff_minmax_rep = { }

		for i in range(len(self.metric_names)):
			all_values_abs = [ residue.metrics[i]['absolute_value'] for residue in self.all_residues if residue.metrics[i]['absolute_value'] is not None ]
			all_values_rep = [ residue.metrics[i]['represented_value'] for residue in self.all_residues if residue.metrics[i]['represented_value'] is not None ]
			if len(all_values_abs) == 0 or len(all_values_rep) == 0:
				raise Exception('One or more metrics have no values specified')
			mean_abs, sdev_abs = sum(all_values_abs)/len(all_values_abs), np.std(all_values_abs)
			mean_rep, sdev_rep = sum(all_values_rep)/len(all_values_rep), np.std(all_values_rep)
			self.metric_stats_abs[i] = (mean_abs, sdev_abs)
			self.metric_stats_rep[i] = (mean_rep, sdev_rep)
			self.metric_quartiles_abs[i] = [ np.percentile(all_values_abs, q) for q in (0, 25, 50, 75, 100) ]
			self.metric_quartiles_rep[i] = [ np.percentile(all_values_rep, q) for q in (0, 25, 50, 75, 100) ]
			metric_sqdiff_polarities_abs = (1 if self.metric_quartiles_abs[i][0] > mean_abs else -1, 1 if self.metric_quartiles_abs[i][-1] > mean_abs else -1)
			metric_sqdiff_polarities_rep = (1 if self.metric_quartiles_rep[i][0] > mean_rep else -1, 1 if self.metric_quartiles_rep[i][-1] > mean_rep else -1)
			self.metric_sqdiff_minmax_abs[i] = (metric_sqdiff_polarities_abs[0]*((self.metric_quartiles_abs[i][0]-mean_abs)/sdev_abs)**2, metric_sqdiff_polarities_abs[1]*((self.metric_quartiles_abs[i][-1]-mean_abs)/sdev_abs)**2)
			self.metric_sqdiff_minmax_rep[i] = (metric_sqdiff_polarities_rep[0]*((self.metric_quartiles_rep[i][0]-mean_rep)/sdev_rep)**2, metric_sqdiff_polarities_rep[1]*((self.metric_quartiles_rep[i][-1]-mean_rep)/sdev_rep)**2)


def _ndist_between(x0, x1, stds=3):
	delta = x1 - x0
	nd = norm.ppf(random(), loc=0.5, scale=0.5/stds)
	nd = min(max(nd, 0), 1)
	# ^ PPF is "percentage-point function", the scipy name for the probit function
	# (aka normal distribution quantile function or inverse CDF)
	y = x0 + nd * delta
	return y


def generate_synthetic(num_chains, chain_length, chain_flaws):
	chains = [ ]
	for i in range(num_chains):
		residues = [ ]
		for j in range(chain_length):
			limits = (0, 30) if j < chain_flaws else (40, 100)
			code = RESIDUES[int(random() * len(RESIDUES))]
			metrics = { }
			metrics.append({ 'name' : 'Metric 1',
							 'shortname' : 'M1',
							 'optimisation' : 'maximise',
							 'absolute_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2),
							 'represented_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2) })
			metrics.append({ 'name' : 'Metric 2',
							 'shortname' : 'M2',
							 'optimisation' : 'maximise',
							 'absolute_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2),
							 'represented_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2) })
			metrics.append({ 'name' : 'Metric 3',
							 'shortname' : 'M3',
							 'optimisation' : 'maximise',
							 'absolute_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2),
							 'represented_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2) })
			metrics.append({ 'name' : 'Metric 4',
							 'shortname' : 'M4',
							 'optimisation' : 'maximise',
							 'absolute_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2),
							 'represented_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2) })
			metrics.append({ 'name' : 'Metric 5',
							 'shortname' : 'M5',
							 'optimisation' : 'maximise',
							 'absolute_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2),
							 'represented_value' : round(_ndist_between(limits[0], limits[1], stds=3), 2) })
			residues.append(ReportResidue(code, metrics))
		shuffle(residues)
		chains.append(ReportChain(i, residues))
	return ReportProtein(chains)


def generate_from_model(metrics_model):
	chains = [ ]
	for i, metrics_chain in enumerate(metrics_model):
		residues = [ ]
		for j, residue in enumerate(metrics_chain):
			metrics = [ ]
			metrics.append({ 'name' : 'B-factor Avg',
							 'shortname' : 'B Av.',
							 'optimisation' : 'minimise',
							 'absolute_value' : round(residue.avg_b_factor, 2),
							 'represented_value' : -round(residue.avg_b_factor, 2) })
			metrics.append({ 'name' : 'B-factor Max',
							 'shortname' : 'B Max',
							 'optimisation' : 'minimise',
							 'absolute_value' : round(residue.max_b_factor, 2),
							 'represented_value' : -round(residue.max_b_factor, 2) })
			metrics.append({ 'name' : 'B-factor Std',
							 'shortname' : 'B Stdev',
							 'optimisation' : 'minimise',
							 'absolute_value' : round(residue.std_b_factor, 2),
							 'represented_value' : -round(residue.std_b_factor, 2) })
			metrics.append({ 'name' : 'Ramachandran',
							 'shortname' : 'Rama Prob.',
							 'optimisation' : 'maximise',
							 'absolute_value' : round(residue.ramachandran_probability, 2) if type(residue.ramachandran_probability)==float else None,
							 'represented_value' : -1/round(residue.ramachandran_probability+0.05, 2) if type(residue.ramachandran_probability)==float else None })
			metrics.append({ 'name' : 'Rotamer',
							 'shortname' : 'Rotamer Prob.',
							 'optimisation' : 'maximise',
							 'absolute_value' : round(residue.rotamer_score, 2) if type(residue.rotamer_score)==float else None,
							 'represented_value' : -round(residue.rotamer_score, 2) if type(residue.rotamer_score)==float else None })
			residues.append(ReportResidue(residue.code, metrics))
		chains.append(ReportChain(i, residues))
	return ReportProtein(chains)


def generate_from_pdb_files(pdb_filenames):
	try:
		from clipper_python import _clipper as clipper
		from clipper_tools.metrics import MetricsModel, MetricsChain, MetricsResidue
	except ImportError:
		raise Exception('Import error. Ensure Clipper-Python and Clipper-tools modules are installed.')
	metrics_models = [ ]
	for filename in pdb_filenames:
		print('Analysing ' + filename)
		fpdb = clipper.MMDBfile()
		fpdb.read_file(filename)
		mmol = clipper.MiniMol()
		fpdb.import_minimol(mmol)
		metrics_model = MetricsModel(mmol)
		metrics_models.append(metrics_model)
	proteins = [ generate_from_model(model) for model in metrics_models ]
	return ReportProteinGroup(proteins)
