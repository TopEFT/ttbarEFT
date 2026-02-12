import os
import argparse
import pprint
import yaml

import pickle
import gzip
import numpy as np
import awkward as ak

import hist
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.utils as utils

import mplhep as hep
import matplotlib.pyplot as plt

def get_sow_SM(h):
    h_EFT = h.as_hist({})
    vals = h_EFT.values()
    return vals

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Customize inputs')
    parser.add_argument('--files', '-f', action='extend', nargs='+', required = True, help = "Specify a list of pkl.gz to run over.")
    args = parser.parse_args()

    # Set variables from command line arguments
    files = args.files

    var_to_save = ['nEvents', 'SumOfWeights', 'SumOfWeights_SM', 'sumw2']
        

    hist_dict = pickle.load(gzip.open(files[0]))['sow']

    saved_vals = {}
    for proc_name in hist_dict.keys():

        hists = hist_dict[proc_name]
        saved_vals[proc_name] = {}

        for var in hists.keys():
            if var in var_to_save:
                vals = get_sow_SM(hists[var])
                saved_vals[proc_name][var] = vals[0]
            else:
                continue

        Neff = np.divide((np.square(saved_vals[proc_name]['SumOfWeights_SM'])), saved_vals[proc_name]['sumw2'])
        saved_vals[proc_name]['Neff'] = Neff

    pprint.pprint(saved_vals)
