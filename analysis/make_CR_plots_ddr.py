import argparse  
import os
import pickle
import gzip
import fnmatch

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
import warnings
warnings.filterwarnings("ignore", message=".*List indexing selection is experimental.*")


plot_title = {"ee_chan": r"$ee$",
              "mm_chan": r"$\mu\mu$",
              "em_chan": r"$e\mu$"}

# label and color based on process axis name, keeps coloring and naming consistent
mc_process_styles = {
    'TTTo2L2Nu*':{'label': r'TTTo2L2Nu', 'color': '#bd1f01'},
    'DY*' :{'label': r'DYJetsToLL', 'color': '#3f90da'},
    'WWTo2L2Nu*':{'label': r'WWTo2L2Nu', 'color': '#b9ac70'},
    'tW*':{'label': r'tW', 'color': '#832db6'},
    'TTGJets*':{'label': r'TTGJets', 'color': '#94a4a2'},
    'ttW*':{'label': r'ttW', 'color': 'pink'},
    'ttZ*':{'label': r'ttZ', 'color': '#a96b59'},
    'WJetsToLNu*':{'label': r'WJetsToLNu', 'color': '#e76300'},
    'WZTo3LNu*':{'label': r'WZTo3LNu', 'color': '#717581'},
    'ZZTo4L*':{'label': r'ZZTo4L', 'color': '#92dadd'},
    'Others': {'label': 'Others', 'color': '#ffa90e'},
}


process_grouping = {
        "TTTo2L2Nu": ["TTTo2L2Nu_centralUL17", "TTTo2L2Nu_centralUL18"],
        "DYJetsToLL": ["DYJetsToLL_centralUL17", "DY10to50_centralUL18", "DY50_centralUL18"],
        "tW": ["tW"]
    }


def get_plotting_style(proc_name, styles_dict=mc_process_styles):
    """
    Checks if proc_name matches a key or a pattern in the styles_dict.
    """
    # Try exact match first
    if proc_name in styles_dict:
        return styles_dict[proc_name]
    
    # Try pattern matching (e.g., 'TTTo2L2Nu*')
    for pattern, style in styles_dict.items():
        if fnmatch.fnmatch(proc_name, pattern):
            return style
            
    # Fallback to making it gray if no color provided
    return {'label': proc_name, 'color': 'gray'}



def group_hist_processes(h, process_map, others_name="Others"):
    """
    Groups a histogram's 'process' axis using a mapping dictionary.
    
    Args:
        h: The original histogram
        process_map: dict of { "NewName": ["OldName1", "OldName2"], ... }
        others_name: Name for any processes not mentioned in the map
    """
    all_procs = list(h.axes['process'])
    
    # Identify which processes from the original hist are vs are NOT in the mapping
    mapped_procs = [p for sublist in process_map.values() for p in sublist]
    remaining_procs = [p for p in all_procs if p not in mapped_procs]

    # Build the list of new categories
    final_categories = list(process_map.keys())
    if remaining_procs:
        final_categories.append(others_name)

    # Initialize the new histogram with the same structure as the old one
    h_grouped = hist.Hist(
        hist.axis.StrCategory(final_categories, name="process", growth=True),
        *[ax for ax in h.axes if ax.name != "process"],
        storage=h.storage_type()
    )

    # Fill the mapped groups
    for new_name, old_names in process_map.items():
        # Only sum processes that actually exist in the current histogram
        existing_olds = [p for p in old_names if p in all_procs]
        if existing_olds:
            # Sum the contributions from the list of old names
            summed_view = h[{"process": existing_olds}][{"process": sum}].view(flow=True)
            h_grouped[{"process": new_name}] = summed_view

    # Group everything else into "Others" process axis
    if remaining_procs:
        others_view = h[{"process": remaining_procs}][{"process": sum}].view(flow=True)
        h_grouped[{"process": others_name}] = others_view

    return h_grouped


def make_cr_fig(h_data, h_mc, var, procs_to_group=None, plot_err=False, h_sumw2=None):

    if procs_to_group:
        h_mc_all = group_hist_processes(h=h_mc, process_map=procs_to_group)
    else: 
        h_mc_all = h_mc

    # create stack for MC plot and make labels/colors based on processes
    mc_stack = h_mc_all.stack("process") #for s in mc_stack: print(s.name)
    mc_labels = [get_plotting_style(s.name, mc_process_styles)['label'] for s in mc_stack]
    mc_colors = [get_plotting_style(s.name, mc_process_styles)['color'] for s in mc_stack]
    # mc_labels = [mc_process_styles.get(s.name, {}).get('label', s.name) for s in mc_stack]
    # mc_colors = [mc_process_styles.get(s.name, {}).get('color', 'gray') for s in mc_stack]

    hep.style.use("CMS")
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1), 'hspace':0.05},
        sharex=True
    )

    ### MC Plot
    hep.histplot(
        list(mc_stack),
        stack=True,
        histtype="fill",
        yerr=False,
        label=mc_labels, #label=[s.name for s in mc_stack],
        color=mc_colors,
        ax=ax
    )

    ### Data Plot
    hep.histplot(
        h_data, 
        histtype='errorbar',
        yerr=True,
        color='black',
        label='Data',
        ax=ax
    )

    ### Ratio Plot 
    h_total_mc = h_mc_all.project(var)
    ratio_hist = h_data / h_total_mc

    # if plot_err:
    #     # ratio_err = (data_vals/mc_vals)*(sqrt((data_uncert/data_vals)^2 + (mc_uncert/mc_vals)^2))
    #     data_vals = h_data.values()
    #     mc_vals = h_total_mc.values()
    #     data_uncert = np.sqrt(h_data.variances())
    #     mc_uncert = np.sqrt(h_total_mc.variances())

    #     # handle zeros
    #     data_term = np.divide(data_uncert, data_vals, out=np.zeros_like(data_vals), where=data_vals>0)
    #     mc_term = np.divide(mc_uncert, mc_vals, out=np.zeros_like(mc_vals), where=mc_vals>0)
    #     ratio_errs = (ratio_hist.values())*np.sqrt(data_term**2 + mc_term**2)
    
    # else: 
    #     ratio_err = False

    hep.histplot(
        ratio_hist,
        yerr=False,
        histtype='errorbar',
        color='black',
        ax=rax,
    )  


    if plot_err:
        ### add mc uncertainty band to ratio plot ###
        mc_vals = h_total_mc.values()
        mc_errs = np.sqrt(h_sumw2.project(f"{var}_sumw2").values())
        bins = h_total_mc.axes[var].edges
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

        # Add to the main plot
        err_m = mc_vals - mc_errs
        err_p = mc_vals + mc_errs

        ax.stairs(
            values=err_m,
            baseline=err_p,
            edges=bins,
            fill=False,
            color='dimgrey',
            # alpha=0.2,
            hatch='\\\\\\\\\\\\\\',
            linewidth= 0,
            label='Stat. Unc.'
        )

        # err_m = np.append(err_m, err_m[-1])
        # err_p = np.append(err_p, err_p[-1])
        # ax.fill_between(bins, err_m, err_p, step='post', facecolor='none', edgecolor='gray', label='Stat err', hatch='\\\\')


    # General formatting
    hep.cms.label("Preliminary", data=True, loc=0, ax=ax)                   #if data=False adds "Simulation" to label
    hep.cms.lumitext(f"(13 TeV)", ax=ax)                                   #TODO add lumi to this
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)   # Scientific notation
    ax.get_yaxis().get_offset_text().set_position((-0.085, 1.05))           # Shift multiplier position out
    # ax.set_ylim(0, ax.get_ylim()[1] * 1.3                                 # or add this, move label inside (loc=2), and put sci ticklabel in default spot

    # Main plot formatting
    ax.set_xlabel("")
    ax.set_ylabel("Events")
    ax.set_xmargin(0)                           # makes 0 on x-axis start at left edge
    if ("eta" in var) or ("phi" in var): 
        leg_loc = 'lower center'
    else: 
        leg_loc = 'upper right'
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc=leg_loc, fontsize=12) # Makes colors on plot appear top to bottom same order as plot 

    # Ratio plot formatting
    rax.set_xlabel(h_total_mc.axes[0].label if h_total_mc.axes[0].label else h_total_mc.axes[0].name)
    rax.set_ylabel("Data/MC")
    rax.set_ylim([0.8, 1.5])                    # rax.set_ylim([0.5, 1.5])
    rax.set_yticks([0.8, 1.0, 1.2, 1.4])        # rax.set_yticks([0.5, 1.0, 1.5]) 
    rax.axhline(y=1.0, color='black', linestyle='--', alpha=0.5)
    rax.set_xmargin(0)                          # makes 0 on x-axis start at left edge
    if "eta" in var: 
        rax.set_xlim([-3, 3])
    elif "phi" in var: 
        rax.set_xlim([-4, 4])


    return fig, ax, rax    


if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description = 'Customize inputs')
    parser.add_argument("--data", required=True, help="path to pkl file of data histograms")
    parser.add_argument("--mc", required=True, help="path to pkl file of MC histograms")
    parser.add_argument("--outdir", default='.', help="output directory")
    parser.add_argument("--outtitle", default='', help="extra title to add to png file names")
    parser.add_argument("--doerror", action='store_true', help="plot uncertainties")
    parser.add_argument("--ylog", action='store_true', help="make plots with log scale on yaxis")

    args = parser.parse_args()
    data_pkl = args.data
    mc_pkl = args.mc 
    outdir = args.outdir
    outtitle = args.outtitle 
    doerror = args.doerror
    ylog = args.ylog

    # make output directory if it doesn't already exist (for running locally)
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    var_to_restrict = ['mll', 'ptll', 'l0pt', 'l1pt']

    hists_data = pickle.load(gzip.open(data_pkl))
    hists_mc = pickle.load(gzip.open(mc_pkl))

    # for channel, ch_dict in hists_data.items():
    for channel, ch_dict in hists_mc.items():

        for var in ch_dict.keys():
            if "sumw2" in var: continue

            h_mc = hists_mc[channel][var].as_hist({})
            h_data = hists_data[channel][var].as_hist({})

            if h_mc.sum() == 0:
                print(f"Skipping {var} in {channel}: MC is empty.")
                continue
            
            if "systematic" in h_mc.axes.name:
                h_mc = h_mc[{"systematic":"nominal"}]

            if "systematic" in h_data.axes.name:
                h_data = h_data[{"systematic":"nominal"}]
            h_data = h_data.project(var)

            # h_mc = hists_mc[channel][var].as_hist({})
            # h_data = hists_data[channel][var].as_hist({}).project(var) # equiv to [{'process':sum}]

            if (h_mc.sum() == 0) or (h_data.sum() == 0):
                print(f"Skipping {var} - No entries found.")
                continue

            if doerror and f"{var}_sumw2" in hists_mc[channel]:
                plot_err = True
                h_sumw2 = hists_mc[channel][f"{var}_sumw2"].as_hist({})
            else:
                plot_err = False 
                h_sumw2 = None

            fig, ax, rax = make_cr_fig(h_data=h_data, h_mc=h_mc, var=var, procs_to_group = process_grouping, plot_err=plot_err, h_sumw2=h_sumw2)
            title = f"{plot_title[channel]}"
            ax.set_title(f"{title}")

            # if var in var_to_restrict:
                # ax.set_xlim([0, 200])
                # rax.set_xlim([0, 200])

            figname = f"{channel}_{var}_CR{outtitle}"
            plt_tools.save_figure(fig, figname, outdir)

            if ylog: 
                ax.set_yscale('log')
                figname = f"{channel}_{var}_CR{outtitle}_log"
                plt_tools.save_figure(fig, figname, outdir)
            
            plt.close(fig)
