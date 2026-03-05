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

def get_binomial_error(h_eff, h_tot):
        """
        Calculates the binomial error
        """
        eff = np.array(h_eff)
        tot = np.array(h_tot)
        
        with np.errstate(divide='ignore', invalid='ignore'):
            error = np.sqrt(eff * (1 - eff) / tot)
            
        return error.tolist()

def make_err_plot(h_val, xbins, ybins, title, hmin=None, hmax=None, ncolors=20):

    # fig, ax = plt.subplots(figsize=(16, 8))
    fig, ax = plt.subplots(figsize=(16, 8))

    if hmin==None:
        hmin = np.nanmin(h_val)
    if hmax==None:
        hmax = np.nanmax(h_val)

    cmap = plt.get_cmap("viridis", ncolors) #18

    hep.hist2dplot(h_val, xbins=xbins, ybins=ybins,labels=True, labels_round=3, labels_fontsize=6, cbarextend=True, cmap=cmap, cmin=hmin, cmax=hmax)
    ax.set_title(f"{title}")
    ax.set_xlim([30, 1000])

    return fig, ax


def make_2d_plot(h_2d, title, hmin=None, hmax=None, ncolors=20):

    # fig, ax = plt.subplots(figsize=(16, 8))
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
    hists = hists_all['jetptetaflav']

    print(hists)
    jetflavors = list(hists.axes['flavour'])

    print(f"jetflavors: {jetflavors}")

    for flav in jetflavors:
        h_num = hists[{'flavour':flav, 'WP':'medium'}]
        h_den = hists[{'flavour':flav, 'WP':'all'}]

        h_eff = h_num / h_den
        # h_eff.values()[:] = np.nan_to_num(h_eff.values(), nan=0.0)

        # hmin = 0.5
        # hmax = 0.88

        fig, ax = make_2d_plot(h_eff, title=f"Eff_{flav}", ncolors=19)
        plt_tools.save_figure(fig, f"Eff_{flav}", outdir)

        fig, ax = make_2d_plot(h_num, title=f"Nbtag_{flav}", ncolors=19)
        plt_tools.save_figure(fig, f"Nbtag_{flav}", outdir)

        fig, ax = make_2d_plot(h_den, title=f"Ntotal_{flav}", ncolors=19)
        plt_tools.save_figure(fig, f"Ntotal_{flav}", outdir)

        pt_bins = [20, 30, 50, 70, 100, 140, 200, 300, 600, 1000]
        eta_bins = [0, 0.6, 1.2, 2.4]

        fig, ax = make_err_plot(get_binomial_error(h_eff.values(), h_den.values()), xbins=pt_bins, ybins=eta_bins, title=f"uncert_{flav}", hmin=0, hmax=0.25)
        plt_tools.save_figure(fig, f"uncert_{flav}", outdir)

        # print(f"h_ref: {h_ref.values()}")
        # print(f"h_denom: {hists[{'process': 'TTTo2L2Nu_centralUL17', 'flav': 'b', 'WP':'all'}].values()}")

        # h_eff_vals = h_ref.values()
        # h_denom_vals = hists[{'process': 'TTTo2L2Nu_centralUL17', 'flav': 'b', 'WP':'all'}].values()

    # # make diff plots
    # h_ref = hists[{'process': 'TTTo2L2Nu_centralUL17', 'flav': 'b', 'WP':'medium'}]/hists[{'process': 'TTTo2L2Nu_centralUL17', 'flav': 'b', 'WP':'all'}]
    # h_sub1 = hists[{'process': 'DYJetsToLL_centralUL17', 'flav':'b', 'WP':'medium'}]/hists[{'process': 'DYJetsToLL_centralUL17', 'flav': 'b', 'WP':'all'}]
    # h_sub2 = hists[{'process': 'tW', 'flav':'b', 'WP':'medium'}]/hists[{'process': 'tW', 'flav': 'b', 'WP':'all'}]

    # fig, ax = make_2d_plot(h_ref-h_sub1, title=f"Powheg bEff - DYJets bEff", hmin=-0.1, hmax=0.12, ncolors=6)
    # plt_tools.save_figure(fig, f"Diff_bEff_DYJets", outdir)
    # fig, ax = make_2d_plot(h_ref-h_sub2, title=f"Powheg bEff - tW bEff", hmin=-0.1, hmax=0.12, ncolors=6)
    # plt_tools.save_figure(fig, f"Diff_bEff_tW", outdir)

