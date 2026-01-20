import functools
import operator 
import gzip 
import pickle
import cloudpickle

from topcoffea.modules.histEFT import HistEFT

orig_hists = pickle.load(gzip.open("MC_CRs_ddr.pkl.gz"))

new_hists = {}

for ch in orig_hists.keys():
	new_hists[ch] = {}

	variables = list(orig_hists[ch][next(iter(orig_hists[ch]))].keys())

	for var in variables: 
		h_list = [orig_hists[ch][dataset][var] for dataset in orig_hists[ch]]

		new_hists[ch][var] = functools.reduce(operator.add, h_list)

with gzip.open("fixed_MC_CRs_ddr.pkl.gz", "wb") as fout:
    cloudpickle.dump(new_hists, fout)

orig_hists = pickle.load(gzip.open("DATA_CRs_ddr.pkl.gz"))

new_hists = {}

for ch in orig_hists.keys():
	new_hists[ch] = {}

	variables = list(orig_hists[ch][next(iter(orig_hists[ch]))].keys())

	for var in variables: 
		h_list = [orig_hists[ch][dataset][var] for dataset in orig_hists[ch]]

		new_hists[ch][var] = functools.reduce(operator.add, h_list)

with gzip.open("fixed_DATA_CRs_ddr.pkl.gz", "wb") as fout:
    cloudpickle.dump(new_hists, fout)