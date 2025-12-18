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


def make_cr_fig(h_data, h_mc, name, title):

    total_mc = h_mc[{'process':sum}]

    centers = h_data.axes.centers[0]
    ratio = np.divide(h_data.values(), total_mc.values())

    # set up figure and colors
    # colors = ['tan', 'tab:gray', 'tab:cyan','tab:pink','tab:green','#7a21dd','#5790fc','#f89c20','#e42536']
    colors = ['#e42536', '#5790fc', 'tab:green', '#f89c20','tab:pink', 'tab:cyan', 'tab:gray', 'tab:brown', 'tan', '#7a21dd']

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
    
    ax.set_title(f"{title}")

    rax.set_ylim([0.6, 1.4])
    rax.set_yticks([0.6, 0.8, 1, 1.2, 1.4])
    
    ax.set_xmargin(0)    # makes 0 on x-axis start at left edge
    rax.set_xmargin(0)   # makes 0 on x-axis start at left edge
    
    if ("eta" in name) or ("phi" in name): 
        ax.legend(loc='lower center', fontsize=14)
    else: 
        ax.legend(loc='upper right', fontsize=14)
    # ax.legend(loc='best', bbox_to_anchor=(1.0, 0.8), fontsize=12)

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

    hists_data = pickle.load(gzip.open(data_pkl))
    hists_mc = pickle.load(gzip.open(mc_pkl))

    channel = "ee_chan"
    var = "ptll"

    mc_process_styles = {
        'TTTo2L2Nu_centralUL17':{'label': r'TTTo2L2Nu', 'color': '#bd1f01'},
        'DYJetsToLL_centralUL17' :{'label': r'DYJetsToLL', 'color': '#3f90da'},
        'WWTo2L2Nu_centralUL17':{'label': r'WWTo2L2Nu', 'color': '#b9ac70'},
        'TW_5f':{'label': r'tW', 'color': '#832db6'},
        'TTGJets_centralUL17':{'label': r'TTGJets', 'color': '#94a4a2'},
        'ttW_centralUL17':{'label': r'ttW', 'color': 'pink'},
        'ttZ_centralUL17':{'label': r'ttZ', 'color': '#a96b59'},
        'WJetsToLNu_centralUL17':{'label': r'WJetsToLNu', 'color': '#e76300'},
        'WZTo3LNu_centralUL17':{'label': r'WZTo3LNu', 'color': '#717581'},
        'ZZTo4L_centralUL17':{'label': r'ZZTo4L', 'color': '#92dadd'},
        'Others': {'label': 'Others', 'color': '#ffa90e'},
    }

    other_processes = ['TTGJets_centralUL17', 'ttW_centralUL17', 'ttZ_centralUL17', 'WJetsToLNu_centralUL17', 'WZTo3LNu_centralUL17', 'ZZTo4L_centralUL17']

    h_mc = hists_mc[channel][var].as_hist({})
    h_data = hists_data[channel][var].as_hist({})
    h_data = h_data.project(var)

    all_mc_procs = list(h_mc.axes['process'])
    sep_mc_procs = [p for p in all_mc_procs if p not in other_processes]
    other_mc_procs = [p for p in all_mc_procs if p in other_processes] 

    h_mc_sep = h_mc[{'process': sep_mc_procs}]  # hist with processes that will be individually plotted

    if other_mc_procs: 
        final_categories = sep_mc_procs + ['Others']

        h_mc_all = hist.Hist(
            hist.axis.StrCategory(final_categories, name= "process", growth=True),
            h_mc.axes[1],
            storage=h_mc.storage_type()
        )

        for proc in sep_mc_procs:
            old_idx = h_mc.axes['process'].index(proc)
            new_idx = h_mc_all.axes['process'].index(proc)
            h_mc_all.view(flow=True)[new_idx] = h_mc.view(flow=True)[old_idx]

        others_hist = h_mc[{'process': other_mc_procs}][{'process':sum}]
        others_idx = h_mc_all.axes['process'].index('Others')
        h_mc_all.view(flow=True)[others_idx] = others_hist.view(flow=True)

    else: 
        h_mc_all = h_mc_sep

    mc_stack = h_mc_all.stack("process") #for s in mc_stack: print(s.name)
    mc_labels = [mc_process_styles.get(s.name, {}).get('label', s.name) for s in mc_stack]
    mc_colors = [mc_process_styles.get(s.name, {}).get('color', 'gray') for s in mc_stack]

    hep.style.use("CMS")
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1), 'hspace':0.05},
        sharex=True
    )

    # MC Plot
    hep.histplot(
        list(mc_stack),
        stack=True,
        histtype="fill",
        label=mc_labels, #label=[s.name for s in mc_stack],
        color=mc_colors,
        ax=ax
    )

    # Data Plot
    hep.histplot(
        h_data, 
        histtype='errorbar',
        color='black',
        label='Data',
        ax=ax
    )

    h_total_mc = h_mc_all.project(var)
    ratio_hist = h_data / h_total_mc
    # # ratio_err = (data_vals/mc_vals)*(sqrt((data_uncert/data_vals)^2 + (mc_uncert/mc_vals)^2))
    data_vals = h_data.values()
    mc_vals = h_total_mc.values()
    data_uncert = np.sqrt(h_data.variances())
    mc_uncert = np.sqrt(h_total_mc.variances())

    # handle zeros
    data_term = np.divide(data_uncert, data_vals, out=np.zeros_like(data_vals), where=data_vals>0)
    mc_term = np.divide(mc_uncert, mc_vals, out=np.zeros_like(mc_vals), where=mc_vals>0)

    ratio_errs = (ratio_hist.values())*np.sqrt(data_term**2 + mc_term**2)

    hep.histplot(
        ratio_hist,
        yerr=None, #can set to ratio_errs
        histtype='errorbar',
        color='black',
        ax=rax,
    )  

    # ### add mc uncertainty band to ratio plot ###
    # mc_vals = h_total_mc.values()
    # mc_errs = np.sqrt(h_total_mc.variances())
    # rel_mc_err = np.divide(mc_errs, mc_vals, out=np.zeros_like(mc_vals), where=mc_vals>0)

    # rax.stairs(
    #     values= 1 + rel_mc_err,      # Upper boundary
    #     baseline= 1 - rel_mc_err,   # Lower boundary
    #     edges=h_total_mc.axes[0].edges,
    #     fill=True,
    #     color='gray',
    #     alpha=0.4,
    #     hatch='////',              # This creates the diagonal lines
    #     label='MC Stat. Unc.'
    # )

    # # Add to the main plot
    # ax.stairs(
    #     values=mc_vals + mc_errs,
    #     baseline=mc_vals - mc_errs,
    #     edges=h_total_mc.axes[0].edges,
    #     fill=True,
    #     color='gray',
    #     alpha=0.4,
    #     hatch='////'
    # )

    # General formatting
    hep.cms.label("Preliminary", data=True, loc=0, ax=ax) #if data=False adds "Simulation" to label
    hep.cms.lumitext(fr"(13 TeV)", ax=ax) #TODO add lumi to this
    ax.ticklabel_format(style="sci", scilimits=(-3, 3), useMathText=True)   # Scientific notation
    ax.get_yaxis().get_offset_text().set_position((-0.085, 1.05))           # Shift multiplier position out
    # ax.set_ylim(0, ax.get_ylim()[1] * 1.3                                 # or add this, move label inside (loc=2), and put sci ticklabel in default spot

    # main plot formatting
    ax.set_ylabel("Events")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='best', fontsize=12) # Makes colors on plot appear top to bottom same order as plot 
    # ax.legend(loc='best', fontsize=12)
    ax.set_xlabel("")

    # ratio plot formatting
    rax.set_xlabel(h_total_mc.axes[0].label if h_total_mc.axes[0].label else h_total_mc.axes[0].name)
    rax.set_ylabel("Data / MC")
    # rax.set_xlabel(var)
    rax.set_ylim([0.5, 1.5])
    # rax.set_yticks([0.5, 1.0, 1.5]) 
    rax.axhline(y=1.0, color='black', linestyle='--', alpha=0.5)


    plt_tools.save_figure(fig, "test", ".")

    # for ch, ch_dict in hists_data.items():  # loop over channels
    #     for dataset, hist_dict in ch_dict.items(): # loop over datasets


    exit()

    plot_title = {"ee_chan": r"$ee$",
                    "mm_chan": r"$\mu\mu$",
                    "em_chan": r"$e\mu$"}

    var_to_restrict = ['mll', 'ptll', 'l0pt', 'l1pt']


    for ch in plot_title.keys(): 
        for name in hists_data: 
            # h_data = hists_data[name][{'process':sum}][{'systematic': "nominal"}].as_hist({})[{'channel':ch}]
            # h_mc = hists_mc[name][{'systematic': "nominal"}].as_hist({})[{'channel':ch}]

            h_data = hists_data[ch][name][{'process':sum}][{'systematic': "nominal"}].as_hist({})
            h_mc = hists_mc[ch][name][{'systematic': "nominal"}].as_hist({})

            fig, ax, rax = make_cr_fig(h_data, h_mc, name, plot_title)

            if name in var_to_restrict:
                ax.sel_xlim([0, 200])
                rax.set_xlim([0, 200])

            figname = ""
            if title is not None: 
                figname += f"{title}_"
            figname += f"{name}_CR_{ch}"

            plt_tools.save_figure(fig, figname, outdir)
            plt.close(fig)
