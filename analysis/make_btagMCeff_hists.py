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


def make_2d_plot(h_2d, title, hmin=None, hmax=None, ncolors=20):

    fig, ax = plt.subplots(figsize=(16, 8))

    if hmin==None:
        hmin = np.nanmin(h_2d.values())
    if hmax==None:
        hmax = np.nanmax(h_2d.values())

    cmap = plt.get_cmap("viridis", ncolors) #18

    hep.hist2dplot(h_2d, ax=ax, labels=True, labels_round=2, labels_fontsize=6, cbarextend=True, cmap=cmap, cmin=hmin, cmax=hmax) #cmin=0, cmax=1.0
    ax.set_title(f"{title}")
    ax.set_xlim([30, 1000])

    return fig, ax 


if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description = 'Customize inputs')
    parser.add_argument("--file", required=True, help='path to pkl file of histograms')
    parser.add_argument("--outdir", default='.', help='path to save figures to')

    args = parser.parse_args()
    hist_file = args.file
    outdir = args.outdir

    # make output directory if it doesn't already exist (for running locally)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    flav_list = ['b', 'c', 'l']

    hists_all = pickle.load(gzip.open(hist_file))['btag']

    # for var, hists in hists_all.items():
    hists = hists_all['jetpteta']
    procs = list(hists.axes['process'])
    jetflavors = list(hists.axes['flav'])

    for p in procs: 
        for flav in jetflavors:
            h_num = hists[{'process':p, 'flav':flav, 'WP':'medium'}]
            h_den = hists[{'process':p, 'flav':flav, 'WP':'all'}]

            h_eff = h_num / h_den
            # h_eff.values()[:] = np.nan_to_num(h_eff.values(), nan=0.0)

            hmin = 0.5
            hmax = 0.88

            fig, ax = make_2d_plot(h_eff, title=f"{p}_Eff_{flav}", hmin=hmin, hmax=hmax, ncolors=19)
            plt_tools.save_figure(fig, f"{p}_Eff_{flav}", outdir)

            fig, ax = make_2d_plot(h_num, title=f"{p}_Nbtag_{flav}", hmin=hmin, hmax=hmax, ncolors=19)
            plt_tools.save_figure(fig, f"{p}_Nbtag_{flav}", outdir)

            fig, ax = make_2d_plot(h_den, title=f"{p}_Ntotal_{flav}", hmin=hmin, hmax=hmax, ncolors=19)
            plt_tools.save_figure(fig, f"{p}_Ntotal_{flav}", outdir)


