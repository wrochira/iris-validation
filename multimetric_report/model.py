from random import random, shuffle

from scipy.stats import norm

from . import METRICS, RESIDUES


class Residue(object):
	def __init__(self, limits=(0, 100)):
		super(Residue, self).__init__()
		self.aa = RESIDUES[int(random() * len(RESIDUES))]
		self.metrics = { }
		for metric in METRICS:
			self.metrics[metric] = round(_ndist_between(limits[0], limits[1], stds=3), 2)


class Chain(object):
	def __init__(self, chain_id, length, num_flaws):
		super(Chain, self).__init__()
		self.chain_id = chain_id
		self.residues = [ Residue(limits=(40, 100)) for _ in range(length-num_flaws) ] + [ Residue(limits=(0, 30)) for _ in range(num_flaws) ]
		shuffle(self.residues)


class Protein(object):
	def __init__(self, chains, chain_length, chain_flaws):
		super(Protein, self).__init__()
		self.chains = [ Chain(chain_id, chain_length, chain_flaws) for chain_id in range(chains) ]


def _ndist_between(x0, x1, stds=3):
	delta = x1 - x0
	nd = norm.ppf(random(), loc=1/2, scale=1/2/stds)
	nd = min(max(nd, 0), 1)
	# ^ PPF is "percentage-point function", the scipy name for the probit function
	# (aka normal distribution quantile function or inverse CDF)
	y = x0 + nd * delta
	return y


def generate(chains, chain_length, chain_flaws):
	return Protein(chains, chain_length, chain_flaws)
