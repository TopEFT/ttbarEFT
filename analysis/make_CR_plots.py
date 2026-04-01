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
        "TTTo2L2Nu": ["TTTo2L2Nu_centralUL16APV", "TTTo2L2Nu_centralUL16", "TTTo2L2Nu_centralUL17", "TTTo2L2Nu_centralUL18"],
        "DYJetsToLL": ["DY10to50_centralUL16APV", "DY50_centralUL16APV", "DY10to50_centralUL16", "DY50_centralUL16", "DYJetsToLL_centralUL17", "DY10to50_centralUL18", "DY50_centralUL18"],
        "tW": ["tW"]
    }

rax_styles = {
    'eta': {'xlim':[-3, 3], 'ylim':[0.6, 1.4], 'yticks':[0.6, 0.8, 1.0, 1.2, 1.4]},
    'phi':  {'xlim':[-4, 4], 'ylim':[0.6, 1.4], 'yticks':[0.6, 0.8, 1.0, 1.2, 1.4]},
    'mll': {'ylim':[0.5, 2]}
}


def get_shape_syst_lst(base_histo, ignore_list=[]):
    # Get unique systematic base names (e.g., "ISR" from "ISRUp")
    all_syst_var_lst = list(base_histo.axes["systematic"])
    syst_var_lst = []
    for name in all_syst_var_lst:
        if name.endswith("Up"):
            base_name = name[:-2] # strip "Up" to get the base systematic name
            if base_name in ignore_list: continue
            if base_name not in syst_var_lst:
                syst_var_lst.append(base_name)
                
    return syst_var_lst


def get_shape_syst_arrs(base_histo, syst_var_lst):
    # Get unique systematic base names (e.g., "ISR" from "ISRUp")
    p_arr_rel_lst = []
    m_arr_rel_lst = []

    for syst_name in syst_var_lst:
        if syst_name == "renormfact": 
            continue

        # Identify relevant samples for this systematic
        h_up = base_histo[{"systematic": syst_name + "Up"}]         # select the "Up" variation syst axis
        relevant_samples_lst = list(h_up.axes["process"])           # get the list of processes that have this syst var

        # Calculate Nominal (n_arr)
        n_arr = base_histo[{"process": relevant_samples_lst, "systematic": "nominal"}]      # select the 'nominal' systematic axis
        n_arr = n_arr[{"process": sum}].values()                                            # sum over the relevant processes

        # handle renorm and fact seperately
        if syst_name in ["renorm", "fact"]:
            continue    # TODO: update for the new hist
            # p_arr_rel, m_arr_rel = get_decorrelated_uncty(syst_name, CR_GRP_MAP, relevant_samples_lst, base_histo, n_arr)
        
        # Calculate Up/Down variations
        else:
            u_arr_sum = base_histo[{"process": relevant_samples_lst, "systematic": syst_name + "Up"}]   # select the samples with an Up variation and the right axis
            u_arr_sum = u_arr_sum[{"process": sum}].values()                                            # sum over all relevant samples 
            
            d_arr_sum = base_histo[{"process": relevant_samples_lst, "systematic": syst_name + "Down"}]
            d_arr_sum = d_arr_sum[{"process": sum}].values()

            # Diff with respect to nominal
            u_arr_rel = u_arr_sum - n_arr
            d_arr_rel = d_arr_sum - n_arr
            
            # Just the ones that increase the yield
            p_arr_rel = np.where(u_arr_rel > 0, u_arr_rel, d_arr_rel)
            m_arr_rel = np.where(u_arr_rel < 0, u_arr_rel, d_arr_rel)

        # Add in quadrature
        p_arr_rel_lst.append(p_arr_rel**2)
        m_arr_rel_lst.append(m_arr_rel**2)

    return [np.sum(m_arr_rel_lst, axis=0), np.sum(p_arr_rel_lst, axis=0)]



def get_proc_plotting_style(proc_name, styles_dict=mc_process_styles):
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
 

def make_cr_fig(h_data, h_mc, var, procs_to_group=None, err_p=None, err_m=None, h_sumw2=None, ylog=False, syst_label="Syst. Unc."):

    plot_syst_err = False
    if (err_p is not None) and (err_m is not None):
        plot_syst_err = True


    if procs_to_group:
        h_mc = group_hist_processes(h=h_mc, process_map=procs_to_group)

    # create stack for MC plot and make labels/colors based on processes
    mc_stack = h_mc.stack("process")                                                            #for s in mc_stack: print(s.name)
    mc_labels = [get_proc_plotting_style(s.name, mc_process_styles)['label'] for s in mc_stack]
    mc_colors = [get_proc_plotting_style(s.name, mc_process_styles)['color'] for s in mc_stack]

    hep.style.use("CMS")

    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize =(10, 10),                                        # figsize=(10,12), # figsize=(11.5,10),
        gridspec_kw={'height_ratios': (4, 1), 'hspace':0.15},     # gridspec_kw={'height_ratios': (3, 1), 'hspace':0.05},
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
        markersize=12,
        yerr=True,
        color='black',
        label='Data',
        ax=ax
    )

    ### Ratio Plot 
    h_total_mc = h_mc.project(var)
    ratio_hist = h_data / h_total_mc

    hep.histplot(
        ratio_hist,
        yerr=False,
        histtype='errorbar',
        markersize=12,
        color='black',
        ax=rax,
    )  

    if plot_syst_err:
        mc_vals = h_total_mc.values()
        bin_edges = h_total_mc.axes[var].edges
        bin_centers = h_total_mc.axes[var].centers

        ### add syst uncertainty band to main plot ###
        p_err_arr = np.where(mc_vals>0,err_p,0)
        m_err_arr = np.where(mc_vals>0,err_m,0)

        ax.fill_between(
            bin_edges, 
            np.append(m_err_arr, m_err_arr[-1]),
            np.append(p_err_arr, p_err_arr[-1]),
            step='post', 
            hatch='\\\\\\\\\\',
            # color='dimgrey',
            alpha=0.1,
            label=syst_label
        )

        ### add syst uncertainty band to ratio plot ###
        p_err_arr_ratio = np.where(mc_vals>0,p_err_arr/mc_vals,1)
        m_err_arr_ratio = np.where(mc_vals>0,m_err_arr/mc_vals,1)

        rax.fill_between(
            bin_edges, 
            np.append(m_err_arr_ratio, m_err_arr_ratio[-1]), 
            np.append(p_err_arr_ratio, p_err_arr_ratio[-1]), 
            step='post', 
            hatch='\\\\\\\\\\',
            # color='dimgrey',
            alpha=0.1, 
        )

    # General formatting
    hep.cms.label("Preliminary", data=True, loc=0, ax=ax)                   #if data=False adds "Simulation" to label
    hep.cms.lumitext(f"(13 TeV)", ax=ax)                                    #TODO add lumi to this
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)   # Scientific notation
    ax.get_yaxis().get_offset_text().set_position((-0.085, 1.05))           # Shift multiplier position out
    # ax.set_ylim(0, ax.get_ylim()[1] * 1.05)                                 # or add this, move label inside (loc=2), and put sci ticklabel in default spot

    if ylog:
        ax.set_yscale('log')
        ax.set_ylim(bottom=0.1, top=ax.get_ylim()[1]*10)
    else: 
        ax.set_ylim(bottom=0, top=ax.get_ylim()[1] * 1.05)

    # Main plot formatting
    ax.set_xlabel("")
    ax.set_ylabel("Events")
    ax.set_xmargin(0)                                                       # makes 0 on x-axis start at left edge
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper right', fontsize=12)  # Makes colors on plot appear top to bottom same order as plot 

    # Ratio plot formatting
    rax.set_xlabel(h_total_mc.axes[0].label if h_total_mc.axes[0].label else h_total_mc.axes[0].name)
    rax.set_ylabel("Data/MC")
    rax.set_xmargin(0)                          # makes 0 on x-axis start at left edge
    rax.axhline(y=1.0, color='black', linestyle='--', alpha=0.5)

    rax.set_ylim([0.5, 1.5])                    # rax.set_ylim([0.5, 1.5])
    rax.set_yticks([0.5, 1.0, 1.5])             # rax.set_yticks([0.8, 1.0, 1.2, 1.4]) 

    for rax_style_key, rax_style in rax_styles.items():
        if rax_style_key in var:
            rax.set(**rax_style)
            break

    return fig, ax, rax   


def make_CR_plots_withsyst(hists_mc, hists_data, var_list=[], syst_list=[], procs_to_group=process_grouping, fig_syst_label='all_syst', ylog=False, do_mcerr=True, plottitle='', figtitle=''):

        if not var_list:                    # if no variables to plot provided, plot all 
            var_list = hists_mc.keys()

        for var in var_list:
            if "sumw2" in var: continue
            
            h_mc = hists_mc[var].as_hist({})
            h_data = hists_data[var].as_hist({})
            
            if h_mc.sum() == 0:
                print(f"Skipping {var}: MC is empty.")
                continue
                
            if "systematic" in h_mc.axes.name:
                h_mc_nom = h_mc[{"systematic":"nominal"}]
            else:
                h_mc_nom = h_mc
            if "systematic" in h_data.axes.name:
                h_data = h_data[{"systematic":"nominal"}]

            h_data = h_data.project(var)
            
            if ('sumw2' in list(h_mc.axes['systematic'])) and do_mcerr: 
                mc_err = h_mc[{'systematic':'sumw2'}][{'process': sum}].values()
                syst_label='Total unc.'
            else: 
                mc_err = 0.0
                syst_label='Syst. unc.'
            
            if not syst_list:               # if no systematics to plot provided, plot all
                syst_list = get_shape_syst_lst(h_mc)

            shape_systs_summed_arr_m , shape_systs_summed_arr_p = get_shape_syst_arrs(h_mc, syst_var_lst=syst_list)
            nom_arr_all = h_mc_nom.project(var).values() 
            p_err_arr = nom_arr_all + np.sqrt(shape_systs_summed_arr_p + mc_err) # This is the upper variation for the main plot 
            m_err_arr = nom_arr_all - np.sqrt(shape_systs_summed_arr_m + mc_err) # the is the lower variation for the main plot

            fig, ax, rax = make_cr_fig(
                h_data=h_data, 
                h_mc=h_mc_nom, 
                var=var, 
                procs_to_group=procs_to_group, 
                err_p=p_err_arr, 
                err_m=m_err_arr, 
                ylog=ylog,
                syst_label=syst_label)

            plt.figtext(0.14, 0.84, fig_syst_label, fontsize=20, fontstyle='italic')  #0.72 
            ax.set_title(plottitle)

            figname = f"{figtitle}_{var}_{fig_syst_label}"
            if ylog: 
                figname += "_log"

            plt_tools.save_figure(fig, figname, outdir)
            plt.close(fig)


if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description = 'Customize inputs')
    parser.add_argument("--data", required=True, help="path to pkl file of data histograms")
    parser.add_argument("--mc", required=True, help="path to pkl file of MC histograms")
    parser.add_argument("--outdir", default='.', help="output directory")
    parser.add_argument("--outtitle", default='', help="extra title to add to png file names")
    parser.add_argument("--doerror", action='store_true', help="plot uncertainties")
    parser.add_argument("--ylog", action='store_true', help="make plots with log scale on yaxis")
    parser.add_argument("--var", )

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

    hist_dict_data = pickle.load(gzip.open(data_pkl))
    hist_dict_MC = pickle.load(gzip.open(mc_pkl))

    for channel, ch_dict in hist_dict_MC.items():
        hists_mc=hist_dict_MC[channel]
        hists_data=hist_dict_data[channel]

        # # plot all variables
        # make_CR_plots_withsyst(
        #     hists_mc=hists_mc, 
        #     hists_data=hists_data, 
        #     var_list=[], 
        #     syst_list=[], 
        #     procs_to_group=process_grouping, 
        #     fig_syst_label='allsyst', 
        #     ylog=False, 
        #     plottitle=plot_title[channel], 
        #     figtitle=channel
        # )

        # make log plots for some variables
        var_list = ['l0pt', 'mll', 'j0pt']
        make_CR_plots_withsyst(hists_mc, hists_data, var_list=var_list, syst_list=[], procs_to_group=process_grouping, fig_syst_label='allsyst', ylog=True, plottitle=plot_title[channel], figtitle=channel)

        # make plots with indiviudal syst variations turned on - does not include mc stat error
        all_systs = get_shape_syst_lst(hists_mc['mll'].as_hist({}))
        var_list = ['njets', 'l0eta', 'mllZ']
        for syst in all_systs: 
            make_CR_plots_withsyst(
                hists_mc=hists_mc, 
                hists_data=hists_data, 
                var_list=var_list, 
                syst_list=[syst], 
                procs_to_group=process_grouping, 
                fig_syst_label=syst, 
                ylog=False, 
                do_mcerr=False,
                plottitle=plot_title[channel], 
                figtitle=channel
            )

        # # make single variable plot
        # var_list = ['njets']
        # make_CR_plots_withsyst(
        #     hists_mc=hists_mc, 
        #     hists_data=hists_data, 
        #     var_list=var_list, 
        #     syst_list=[], 
        #     procs_to_group=process_grouping, 
        #     fig_syst_label='allsyst', 
        #     ylog=False, 
        #     plottitle=plot_title[channel], 
        #     figtitle=f"test_{channel}"
        # )


