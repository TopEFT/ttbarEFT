import os

import numpy as np

import hist
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.utils as utils

import mplhep as hep
import matplotlib.pyplot as plt


def calc_sow_array(h, wc_range, wc_name, st_pt=None):
    '''
    Calculate the sum of event weights for a range of values of a 
    single wc while all others are frozen at another value. 
    
    Parameters
    ----------
    h : single histEFT histogram (based on scikithep hist)
        needs to be a 1 bin histogram in order to get the correct output
    wc_range : list
        List of wc values to calculate histogram event weight at
        eg. np.arange(min, max, increment)
    wc_lst : list
        list of all WC contained in the histogram
    wc_name : str
        single wc that will be scanned (the rest are set to st_pt)
    st_pt: dict
        dictionary of wc values to set the rest to

    Returns: 
        weight_pts: list
        [[wc-value], [hist.values() at each wc value]]
    '''

    wc_lst = h._wc_names.keys()

    weight_pts = [] 
    weights = []
    wc_vals = {}

    # fill wc_vals dict
    if st_pt is None: 
        # if no st_pt provided, assume SM (all WCs=0)
        for item in wc_lst: 
            wc_vals[item]=0.0
    else: 
        # check that the st_pt is a subste of wc_lst: 
        assert set(st_pt.keys()).issubset(set(wc_lst)), "Provided dict of wc values should be a subset of the wc_lst in the histEFT."
        #fill wc_vals with values from st_pt. If there is no provided value, set it to 0. 
        for name in wc_lst.keys():
            if name in st_pt.keys():
                wc_vals[name] = st_pt[name]
            else: 
                wc_vals[name] = 0.0

    # loop through different values of the wc coeff, get sum of event weights
    for i in wc_range:
        wc_vals.update({wc_name:i})
        weight = h.as_hist(wc_vals).values()
        assert len(weight)==1, f"Histogram should only contain one bin. len(hist.values)={len(weight)}"
        weights.extend(weight[0])

    weight_pts.append(wc_range)
    weight_pts.append(weights)
        
    # print("in calc_sow_array: \n")
    # print(f"wc_name = {wc_name}")
    # print(f"weight_pts = {weight_pts}")

    return weight_pts

def get_single_hist(fname, hist_name):
    '''Get single histogram from pkl file based on hist name'''
    hists = utils.get_hist_from_pkl(fname, allow_empty=False)
    h = hists[hist_name]
    return h

def make_1d_quad_plot_with_scatter(samples, scatter_lst, wc_max, wc_name, ylabel=r"$\sigma_{SMEFT} /\ \sigma_{SM}$"):
    '''
    Create and save 1d quadratic plot of wc value versus total event weight

    Parameters
    ----------
    samples : dict
        dictionary containing sample names and the weights array (output from calc_sow_array)
        {sample name: [[wc-value],[xsec calculated by histEFT]]}
    figname : str
        Figure name used to save plot (including path if needed)
    wc_max : int
        max wc value scanned, used to set the x-axis range
    wc_name : str
        single wc that is being scanned, used to set x-axis label
    scatter_lst : list 
        list of wc values, xsec values, and xsec uncertainties calculated by MG standalone 
        ouptut of read_MGstandalone_txt() for the given wc being scanned
        [[x-values],[y-values],[uncertainty on y-values]]

    Returns 
    -------
    fig, ax : matplotlib figure and axes objects
    '''

    hep.style.use("CMS")
    params = {'axes.labelsize': 20,
          'axes.titlesize': 20,
          'legend.fontsize':20}
    plt.rcParams.update(params)
    fig, ax = plt.subplots()

    ax.scatter(scatter_lst[0], scatter_lst[1], label = "MG standalone")
    ax.errorbar(scatter_lst[0], scatter_lst[1], yerr = scatter_lst[2], xerr = None, capsize=5, ls='none')

    for item in samples:
        ax.plot(samples[item][0], samples[item][1], label=str(item))
    
    # ax.legend(loc='upper right', fontsize='medium')
    ax.legend(loc=(1.04, 0.5), fontsize='medium')
    ax.set_xlim([-wc_max, wc_max])
    ax.set_xlabel(wc_name, fontsize = 'large')
    ax.set_ylabel(ylabel, fontsize = 'large')

    return fig, ax

def read_MGstandalone_txt(fname):
    '''
    Return a dictionary of xsec values at different WC points. 

    Parameters
    ---------- 
    fname : str
        file name (and if not in the same directory, the path) of the txt file that contains 
        wc values and the cross section calculated by stand alone MG. Formatted like
        "ctGIm=-0.7
        5.703 +- 0.01831 pb"

    Returns
    -------
    dedicated: dict
        dictionary of wc names and xsec points with uncertainties
        {wc name: [[wc_values], [xsec_values], [xsec_uncertainties]]}
    '''
    
    with open(fname, "r") as f: 
        lines = f.read().splitlines() 
    f.close() 

    dedicated = {}
    
    i = 0
    while i < len(lines):
        first = lines[i].split("=")
        second = lines[i+1].split("+-")

        name = first[0]
        xval = float([first[1]][0])
        yval = float([second[0]][0])
        ysigma = float(second[1][:-3])

        if name not in dedicated.keys():
            dedicated[name] = [[xval], [yval], [ysigma]]
        else: 
            dedicated[name][0].append(xval)
            dedicated[name][1].append(yval)
            dedicated[name][2].append(ysigma)   
        i += 2

    return dedicated


def save_figure(fig, figname, outdir=''):
    outname=os.path.join(outdir, f"{figname}.png")
    fig.savefig(outname, bbox_inches='tight')
    print(f'plot saved to {outname}')
    plt.close(fig)
