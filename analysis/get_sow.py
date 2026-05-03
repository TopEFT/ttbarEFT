import os
import argparse

# import pickle
# import gzip
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
        
    for file in files: 
        print(f"FILENAME: {file} \n")
        hists = utils.get_hist_from_pkl(file, allow_empty=False)
        for name in hists.keys():
            h_main=hists[name]
            axes=h_main.axes[0]
            print(f"HISTOGRAM: {name}")
            for ax in axes: 
                h_temp = h_main[ax]
                vals = get_sow_SM(h_temp)
                print(f"    {ax}={vals}")
            print("\n")

