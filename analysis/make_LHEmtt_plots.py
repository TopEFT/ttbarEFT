import argparse  
import os
import pickle
import gzip

import awkward as ak
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from cycler import cycler
import mplhep as hep

import hist
from hist import Hist

from topcoffea.modules.histEFT import HistEFT
from topcoffea.modules import utils
import topcoffea.modules.eft_helper as efth

from ttbarEFT.modules import plotting_tools_histEFT as plt_tools

np.seterr(divide='ignore', invalid='ignore')


def make_plot(hists, procs):

    hep.style.use("CMS")
    fig, ax = plt.subplots(figsize=(12, 12))

    for p in procs: 
        hep.histplot(
            hists[{'process':p}], 
            histtype="fill", 
            yerr=False,
            label=str(p),
            ax=ax
        )

    ax.set_ylabel("Events")
    ax.set_xlabel("mtt")
    ax.legend(loc='best', fontsize=14)
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)   # Scientific notation

    return fig, ax



if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description = 'Customize inputs')
    parser.add_argument("--hists", required=True, help="path to pkl file of data histograms")
    parser.add_argument("--outdir", required=False, default='./', help='path for output directory')

    args = parser.parse_args()
    hist_file = args.hists
    outdir = args.outdir

    # make output directory if it doesn't already exist (for running locally)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    hists= pickle.load(gzip.open(hist_file))

    proc_list = list(hists['mtt']['mtt'].axes['process'])
    years = ['UL16_', 'UL16APV_', 'UL17', 'UL18']

    wcs = {
      'SM': {}, 
      'StartingPoint': {'cQd1': 1.5, 'ctj1': 1.5, 'cQj31': 1.5, 'ctj8': 1.5, 
                'ctd1': 1.5, 'ctd8': 1.5, 'ctGRe': -0.5, 'ctGIm': -0.5, 
                'cQj11': 1.5, 'cQj18': 1.5, 'ctu8': 1.5, 'cQd8': 1.5, 
                'ctu1': 1.5, 'cQu1': 1.5, 'cQj38': 1.5, 'cQu8': 1.5},
      'LargeWC': {'cQd1': 3, 'ctj1': 3, 'cQj31': 4, 'ctj8': 4, 
                'ctd1': 3, 'ctd8': 3, 'ctGRe': 1, 'ctGIm': 1, 
                'cQj11': 2, 'cQj18': 4, 'ctu8': 4.5, 'cQd8': 4.5, 
                'ctu1': 2, 'cQu1': 2, 'cQj38': 3.5, 'cQu8': 2.5},
      'VeryLargeWC': {'cQd1': 8, 'ctj1': 10, 'cQj31': 12, 'ctj8': 11, 
                    'ctd1': 9, 'ctd8': 9, 'ctGRe': -5, 'ctGIm': 10, 
                    'cQj11': 10, 'cQj18': 11, 'ctu8': 9, 'cQd8': 13, 
                    'ctu1': 12, 'cQu1': 13, 'cQj38':14, 'cQu8': 12}
    }

    for year in years: 
        year_procs = [p for p in proc_list if year in p]

        for rwgt in wcs.keys():
        # for rwgt in ['VeryLargeWC']:
            fig, ax = make_plot(hists=hists['mtt']['mtt'].as_hist(wcs[rwgt]), procs=year_procs)
            ax.set_title(f"Reweighted to {rwgt}")
            figname = f"mtt_{rwgt}_{year}"
            plt_tools.save_figure(fig, figname, outdir)

