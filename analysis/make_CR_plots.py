import argparse  
import os

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


def make_cr_fig(h_data, h_mc):

    total_mc = h_mc[{'process':sum}]

    centers = h_data.axes.centers[0]
    ratio = np.divide(h_data.values(), total_mc.values())

    # set up figure and colors
    colors = ['tan', 'tab:gray', 'tab:cyan','tab:pink','tab:green','#7a21dd','#5790fc','#f89c20','#e42536']

    hep.style.use("CMS")
    # Initialize figure and axes
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)
    ax.set_prop_cycle(cycler(color=colors)) 

    # plot histograms
    h_data.plot1d(ax=ax, yerr=False, histtype='errorbar', linewidth=2, label='Data', color='black')
    h_mc.plot1d(ax=ax, yerr=False, stack=True, histtype='fill', linewidth=2)

    # plot ratio 
    rax.scatter(centers, ratio, color='black')

    # formatting
    xlabel = h_data.axes.label[0]
    ax.set_ylabel("Events")
    ax.set_xlabel("")
    rax.set_xlabel(xlabel)
    rax.set_ylabel("Data/MC")
    
    ax.set_title(f"{ch}")
    ax.legend(loc='best', ncol=2, fontsize=12)

    rax.set_ylim([0.6, 1.4])
    rax.set_yticks([0.6, 0.8, 1, 1.2, 1.4])
    
    ax.set_xmargin(0)    # makes 0 on x-axis start at left edge
    rax.set_xmargin(0)   # makes 0 on x-axis start at left edge
    
    if "eta" in name: 
        rax.set_xlim([-3, 3])
    elif "phi" in name: 
        rax.set_xlim([-4, 4])
        
    rax.axhline(y=1.0, color='gray', linestyle='--')
    rax.axhline(y=1.2, color='gray', linestyle='--')
    rax.axhline(y=0.8, color='gray', linestyle='--')
    rax.axhline(y=0.6, color='gray', linestyle='--')

    return fig, ax, rax


if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description = 'Customize inputs')
    parser.add_argument("--data", required=True, help="path to pkl file of data histograms")
    parser.add_argument("--mc", required=True, help="path to pkl file of MC histograms")
    parser.add_argument("--outdir", default='.', help="output directory")
    parser.add_argument("--title", default=None, help="extra title to add to png file names")

    args = parser.parse_args()
    data_pkl = args.data
    mc_pkl = args.mc 
    outdir = args.outdir
    title = args.title 

    # make output directory if it doesn't already exist (for running locally)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    hists_data = utils.get_hist_from_pkl(data_pkl, allow_empty=False)
    hists_mc = utils.get_hist_from_pkl(mc_pkl, allow_empty=False)

    channels = ["ee", "mm"]

    for ch in channels: 
        for name in hists_data: 
            h_data = hists_data[name][{'process':sum}].as_hist({})[{'channel':ch}]
            h_mc = hists_mc[name].as_hist({})[{'channel':ch}]

            fig, ax, rax = make_cr_fig(h_data, h_mc)

            figname = ""
            if title is not None: 
                figname += f"{title}_"
            figname += f"{name}_CR_{ch}"

            plt_tools.save_figure(fig, figname, outdir)
