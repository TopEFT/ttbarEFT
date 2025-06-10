import numpy as np
import awkward as ak
import uproot
import os

import hist

import mplhep as hep
import matplotlib.pyplot as plt
import matplotlib as mpl

from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.utils as utils
import ttbarEFT.modules.plotting_tools_histEFT as plotTools


def get_ratio_uncertainty(num_hist, denom_hist):
    '''
    Calculates the propagated uncertainty per bin on the ratio of two historgams

    Parameters
    ----------
        num_hist (scikithep hist): numerator histogram
        denom_hist (scikithep hist): Denominator histogram
    
    Returns:
        list of uncertainties, one entry per histogram bin
    '''

    xvals = num_hist.axes.centers[0]
    yvals_num = num_hist.values()
    yvals_denom = denom_hist.values()
    sigma_num = np.sqrt(num_hist.variances())
    sigma_denom = np.sqrt(denom_hist.variances())

    ratio = np.divide(yvals_num, yvals_denom)

    # calculation for error propagation for ratio = yavls_num/yvals_denom
    # generally, z=x/y; sigma_z = abs(z)sqrt((sigma_x/x)^2+(sigma_y/y)^2)
    sigma_y = np.multiply(np.abs(ratio), np.sqrt(np.add(np.square(np.divide(sigma_num, yvals_num)), np.square(np.divide(sigma_denom, yvals_denom)))))

    return sigma_y


def get_ratio_points(num_hist, denom_hist):
    '''
    Calculates the ratio between two histograms

    Parameters
    ----------
        num_hist (scikithep hist): numerator histogram
        denom_hist (scikithep hist): Denominator histogram

    Returns:
        centers: list of the x-axis point that the ratio corresponds to
        ratio: list of the ratio value (one for each hist bin)
    '''

    num = num_hist.values()
    centers = num_hist.axes.centers[0]
    denom = denom_hist.values()
    ratio = np.divide(num, denom)

    return centers, ratio


def get_rwgt_xsec_list(file, wc_name, wc_max=6.0, hist_name='sow_norm'):
    '''
    Creates a nested list that contains the x and y values to use in a cross section plot
    This is done by reweighting a histEFT histogram to a variety of values for the same WC

    Parameters
    ----------
    file(str): absolute path to historgam file
    wc_name(str): wc name that's being scanned
    hist_name(str): name of the histogram that will be reweighted (e.g. 'sow', 'sow_norm', 'nEvents')
    wc_max(int): maximum value to calculate of the WC

    Returns: weights(list)
        [[wc-value], [hist.values() at each wc value]]
    '''
    
    samples = {}

    if wc_name == 'ctGRe' or wc_name == 'ctGIm':
        wc_range = np.arange(-1.5, 1.5, 0.2)
    else:
        wc_range = np.arange(-wc_max, wc_max+0.5, 0.5)
    h = plotTools.get_single_hist(file, hist_name)
    norm = h.as_hist({}).values()[0] #get SM xsec of the sample and use this for normalization
    weights = plotTools.calc_sow_array(h, wc_range, wc_name)

    if norm != 1.0:
        weights[1] = np.divide(weights[1], norm)

    return weights


def get_scatter_xsec_list(info, norm, norm_uncert):
    '''
    Creates a nested list that contains the x and y values for a plot and the uncertainty on the y value
    This formats the data from a MG standalone xsec calculation into lists that can be plotted

    Parameters
    ----------
    info (nested list): list of the format [[WC values], [xsec corresponding to the wc vals], [uncertainty for each xsec]]
        this is the output from plotTools.read_MGstandalone_txt
    norm (float): this is the number that the xsec is being normalized to (usually the SM cross section)
    norm_uncert(float): this is the uncertainty on the normalization number

    Returns: nested list 
        [x-axis values, y-axis values (the normalized cross section), uncertainties on the normalized cross section]
    '''
    scatter_xvals = info[0]
    scatter_yvals = np.divide(np.array(info[1]), norm)
    scatter_sigma = np.array(info[2])
    sigma_y = np.multiply(scatter_yvals, (np.sqrt(np.add(np.square(np.divide(scatter_sigma, info[1])),np.square(np.divide(norm_uncert, norm))))))

    return [scatter_xvals, scatter_yvals, sigma_y]


def make_1D_quad_plot(samples, wc_name, const=1.0):
    '''
    Create 1D quadratic EFT parameterization plots

    Parameters:
    -----------
        samples (dict): dictionary of samples to plot 
            each sample is contained in a dictionary entry of the following format: 
            {'sample_name': {'data': sample_list, 'type': '', 'label': 'plot label here'}}
                'sample_dict'=the output of either get_rwgt_xsec_list or get_scatter_xsec_list
                'type'='plot' or 'scatter' depending on the contents of sample_dict. This will determine if `ax.plot` or `ax.scatter` is used
    Returns: 
        figure, ax of the 1d quad EFT plot 
    '''

    hep.style.use("CMS")
    fig, ax = plt.subplots()

    for item in samples:
        if samples[item]['type']=='scatter':
            ax.scatter(samples[item]['data'][0], samples[item]['data'][1]*const, label = samples[item]['label'])            
            ax.errorbar(samples[item]['data'][0], samples[item]['data'][1]*const, yerr = samples[item]['data'][2], xerr = None, capsize=5, ls='none')
        elif samples[item]['type']=='plot':
            ax.plot(samples[item]['data'][0], samples[item]['data'][1], label=samples[item]['label'])

    ax.legend(loc='best', fontsize='medium')
    ax.set_xlabel(wc_name, fontsize = 'large')
    ax.set_ylabel(r"$\sigma_{SMEFT} /\ \sigma_{SM}$", fontsize='large')

    return fig, ax 


def make_histEFT_ratio_plot_multi(h_dict, h_denom, ratio_points=True):
    '''
    Parameters: 
        h_dict (dict): dictionary of samples to be plotted. 
                Each will be plotted and then the ratio between it and h_denom will be plotted
                {'hist_name': {'h':histEFT that has been reweighted as histEFT.as_hist({}), 
                                'plot_options': {plot options dict including 'label'}}, ...}
        h_denom (dict): same format as a single entry in h_dict
                this will be used as the denominator for the ratio plot
                {'h':histEFT that has been reweighted as histEFT.as_hist({}), 
                    'plot_ops': {plot options dict including 'label'}}
    '''

    hep.style.use("CMS")
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)
    
    bin_widths = h_denom['h'].axes[1].edges
    centers = h_denom['h'].axes[1].centers
    val_denom = h_denom['h'].values()[0]

    # loop through provided reweight points and plot the historgam and ratio for each
    for name, hist_info in h_dict.items(): 
        h = hist_info['h']
        ratio = h.values()[0]/val_denom

        h.plot1d(ax=ax, **hist_info['plot_options'])
        if ratio_points:
            rax.scatter(centers, ratio, label=hist_info['plot_options']['label'])
        else:
            hep.histplot(ratio, bin_widths, ax=rax)

    # plot the denominator histogram
    h_denom['h'].plot1d(ax=ax, color='black', **h_denom['plot_options'])

    ax.legend(loc='upper right')
    rax.axhline(y=1.0, color='gray', linestyle='--')

    return fig, ax, rax


def make_histEFT_plot_multi(h_dict, **kwargs):
    '''
    Make a figure that contains multiple histograms reweighted to the same point plotted together

    Parameters: 
        h_dict (dict): {'hist_name': {'h':single_histEFT, 
                                      'rwgt': {reweight dictionary}, 
                                      'plot_ops': {plot options dict including 'label'}}, ...}
        **kwargs: keyword arguments. All will be passed to ax.set()
    '''

    hep.style.use("CMS")
    # hep.cms.label(label='Work in progress', com=13)
    fig, ax = plt.subplots()

    # default options that will be overridden by function input if they are given
    default_ax_options = {
        'ylabel': 'Events', 
    }
    ax_options = default_ax_options.copy()
    ax_options.update(**kwargs)

    # loop through the historgams provided and plot them on the same axis
    for name, hist_info in h_dict.items(): 
        h = hist_info['h'].as_hist(hist_info['rwgt'])
        h.plot1d(ax=ax, **hist_info['plot_options'])

    # formatting 
    ax.set(**ax_options)
    ax.legend(loc='upper right')

    return fig, ax


def make_hist_plot_single(h, **kwargs):
    '''
    Returns a figure that contains a single histogram

    Parameters: 
        h: a scikithep histogram
        **kwargs: keyword arguments. The function expects 
            `plot_options` (dict): kwargs passed to hep.histplot() -e.g. stack=False, yerr=True
            `ax_options` (dict): kwargs passed to ax.set() - e.g. xlim, ylim, title, xlabel, ylabel
    '''
    input_plot_options = kwargs.pop('plot_options', {})
    input_ax_options = kwargs.pop('ax_options', {})

    hep.style.use("CMS")
    # hep.cms.label(label='Work in progress', com=13)
    fig, ax = plt.subplots()

    # default options that will be overridden by function input if they are given
    default_ax_options = {
        'ylabel': 'Events', 
    }
    default_plot_options = {
        'stack': False,
        'yerr' : False, 
        'linewidth': 2
    }
    plot_options = default_plot_options.copy()
    plot_options.update(input_plot_options)
    ax_options = default_ax_options.copy()
    ax_options.update(input_ax_options)

    # fill plot 
    hep.histplot(h, ax=ax, **plot_options)

    # formatting 
    ax.set(**ax_options)
    ax.legend(loc='upper right')

    return fig, ax


def make_hist_plot_multi(h_dict, **kwargs):
    '''
    Make a figure that contains multiple histograms plotted together

    Parameters: 
        h_dict (dict): {'hist_name': {'h':single_hist, 'plot_ops': {plot options dictionary}}, ...}
        **kwargs: keyword arguments. Everything will be used for axis options using ax.set()
    '''

    hep.style.use("CMS")
    # hep.cms.label(label='Work in progress', com=13)
    fig, ax = plt.subplots()

    # default options that will be overridden by function input if they are given
    default_ax_options = {
        'ylabel': 'Events', 
    }
    ax_options = default_ax_options.copy()
    ax_options.update(**kwargs)

    # loop through the historgams provided and plot them on the same axis
    for hist in h_dict: 
        hep.histplot(h_dict[hist]['h'], ax=ax, **h_dict[hist]['plot_options'])

    # formatting 
    ax.set(**ax_options)
    ax.legend(loc='upper right')

    return fig, ax


def make_hist_ratio_plot_single(h_num, h_denom, label_num, label_denom, xlabel, title, ratio_hlines=[1,5], ratio_ylim=[0,2]):

    # get ratios and uncertainties
    centers, ratio = get_ratio_points(h_num, h_denom)
    uncert = get_ratio_uncertainty(h_num, h_denom)

    hep.style.use("CMS")
    # hep.cms.label(label='Work in progress', com=13)

    # Initialize figure and axes
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)

    # Plot histograms with errors
    hep.histplot(h_denom, ax=ax, stack=False, yerr=True, linewidth=2, label=label_denom)
    hep.histplot(h_num, ax=ax, stack=False, yerr=True, linewidth=2, label=label_num)
    # Plot ratio with errors
    rax.scatter(centers, ratio)
    rax.errorbar(centers, ratio, xerr = None, yerr = uncert, capsize=5, ls='none')

    ## Formatting
    ax.legend(loc = 'upper right')
    ax.set_ylabel("Events", fontsize='medium')
    ax.set_xlabel("")
    rax.set_ylabel("Ratio", fontsize='medium')
    rax.set_xlabel(xlabel, fontsize="medium")
    rax.set_ylim(ratio_ylim)
    for line in ratio_hlines:
        rax.axhline(y=line, color='gray', linestyle='--')
    plt.figtext(0.13, 0.9, title, fontsize='medium')

    return fig, ax, rax


def save_figure(fig, figname, outdir=''):
    outname=os.path.join(outdir, f"{figname}.png")
    fig.savefig(outname, bbox_inches='tight')
    print(f'plot saved to {outname}')
    plt.close(fig)


def main():
    savedir = '/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tW_plots/fullchecks/'

    ######################
    # Event Weight Plots #
    ######################
    tWLO_stpt2_weights = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_weights.pkl.gz", allow_empty=False)
    
    event_weight_ax_ops = {'xlabel': 'log(event weight)', 'ylabel': 'Counts', 'xlim': [-2.5, 0.1], 'title': 'tWLO_stpt2'}
    event_weight_hist_dict = {}
    for h in tWLO_stpt2_weights: 
        event_weight_hist_dict[h] = {'h': tWLO_stpt2_weights[h], 'plot_options':{'stack':False, 'yerr':True, 'linewidth':2, 'label': h[:-4][8:]}}
    fig, ax = make_hist_plot_multi(event_weight_hist_dict, **event_weight_ax_ops)
    save_figure(fig, "event_weight_stpt2", outdir=savedir)

    # # ### Individual Event Weight Plots ### 
    # # # for h in tWLO_stpt2_weights: 
    # # #     hist_label = h[:-4][8:]
    # # #     ax_ops = {'xlabel': 'log(event weight)', 'ylabel': 'Counts', 'xlim': [-2.5, 0.1]}
    # # #     plot_ops = {'stack':False, 'yerr':True, 'linewidth':2, 'label': hist_label}
    # # #     fig, ax = make_hist_plot_single(tWLO_stpt2_weights[h], ax_options=ax_ops, plot_options=plot_ops)


    ########################################
    # Event Weight Plots with LHE mmnl>100 #
    ########################################

    tWLO_stpt2_weights_LHEmmnl = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_weights_mmnl.pkl.gz", allow_empty=False)
    event_weight_ax_ops = {'xlabel': 'log(event weight)', 'ylabel': 'Counts', 'xlim': [-2.5, 1.0], 'title': 'LHE mmnl>100'}
    event_weight_hist_dict = {}
    for h in tWLO_stpt2_weights_LHEmmnl: 
        event_weight_hist_dict[h] = {'h': tWLO_stpt2_weights_LHEmmnl[h], 'plot_options':{'stack':False, 'yerr':True, 'linewidth':2, 'label': h[:-4][8:]}}
    fig, ax = make_hist_plot_multi(event_weight_hist_dict, **event_weight_ax_ops)
    save_figure(fig, "event_weights_mmnl100", outdir=savedir)


    ######################################
    # hist plots with uncertainties #
    ######################################

    # comparing EFT to powheg at SM 
    tWtop_hist_norm1 = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWtop_powheg_kin_norm.pkl.gz", allow_empty=False)
    tWantitop_hist_norm1 = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWantitop_powheg_kin_norm.pkl.gz", allow_empty=False) 
    tWLO_stpt2_hist_norm1 = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_kin_norm.pkl.gz", allow_empty=False) 

    for name in tWLO_stpt2_hist_norm1: 
        hpowheg = tWtop_hist_norm1[name]+tWantitop_hist_norm1[name]
        hEFT = tWLO_stpt2_hist_norm1[name]
        fig, ax, rax = make_hist_ratio_plot_single(hEFT, hpowheg, "tW_EFT", "tW_powheg", name, "Reweighted to SM, normalized to 1.0")
        # at the SM, the mt2 bins after 125 are basically zero so just cut these out for the powheg comparison
        if name =='mt2':
            ax.set_xlim([0, 125])
        save_figure(fig, f"{name}_norm1", outdir=savedir)


    ######################################
    # histEFT: reweighted mT2 with ratio #
    ######################################

    tWLO_stpt1_histEFT_taus = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt1_histEFT_taus.pkl.gz", allow_empty=False)
    tWLO_stpt2_histEFT_taus = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_histEFT_taus.pkl.gz", allow_empty=False)

    h = tWLO_stpt2_histEFT_taus['mt2']
    h_SM = h.as_hist({})
    rwgt1_dict = {'cQl311=0.5':{'cQl311':0.5}, 'cleQt3Re11=0.5':{'cleQt3Re11':0.5}, 'cleQt1Re11=1.0':{'cleQt1Re11':1.0}}
    rwgt2_dict = {'cQl322=0.5':{'cQl322':0.5}, 'cleQt3Re22=0.5':{'cleQt3Re22':0.5}, 'cleQt1Re22=1.0':{'cleQt1Re22':1.0}}
    rwgt3_dict = {'cQl333=0.5':{'cQl333':0.5}, 'cleQt3Re33=0.5':{'cleQt3Re33':0.5}, 'cleQt1Re33=1.0':{'cleQt1Re33':1.0}}

    plots_to_make = {'1st generation lepton operators': rwgt1_dict, '2nd generation lepton operators': rwgt2_dict, '3rd generation lepton operators': rwgt3_dict}
    
    for name, rwgt_dict in plots_to_make.items():

        plot_options = {'histtype': 'step', 'stack':False, 'yerr':False, 'linewidth':2, 'label':''}
        h_denom = {'h':h_SM, 'plot_options': {**plot_options, 'label':'SM'}}
        h_dict={}

        for rwgt_name, rwgt_val in rwgt_dict.items():
            h_dict[rwgt_name]={'h':h.as_hist(rwgt_val), 'plot_options':{**plot_options, 'label':rwgt_name}}

        fig, ax, rax = make_histEFT_ratio_plot_multi(h_dict, h_denom, ratio_points=False)

        ax_options = {'xlabel':'', 'ylabel':'Events/Bin', 'yscale':'log', 'title':name}
        rax_options = {'xlabel': 'mT2 [GeV]', 'ylabel': 'EFT/SM', 'ylim': [0,13]}
        rax.axhline(y=5.0, color='gray', linestyle='--')
        ax.set(**ax_options)
        rax.set(**rax_options)
        
        save_figure(fig, f"mt2_rwgt_{name[:3]}", outdir=savedir)


    #####################################################################
    # EFT Hists reweighted to different starting point with uncertainty #
    #####################################################################

    tWLO_stpt1_hist_rwgt1 = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt1_hist_rwgt1.pkl.gz", allow_empty=False)
    tWLO_stpt2_hist_rwgt1 = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_hist_rwgt1.pkl.gz", allow_empty=False)
    tWLO_stpt1_hist_rwgt2 = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt1_hist_rwgt2.pkl.gz", allow_empty=False)
    tWLO_stpt2_hist_rwgt2 = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_hist_rwgt2.pkl.gz", allow_empty=False) 

    for name in tWLO_stpt1_hist_rwgt1: 
        h1 = tWLO_stpt1_hist_rwgt1[name]
        h2 = tWLO_stpt2_hist_rwgt1[name]

        fig, ax, rax = make_hist_ratio_plot_single(h2, h1, "Sample2", "Sample1", name, "Reweighted to Starting Point of Sample1")
        rax.set_ylabel('Ratio(S2/S1)')
        save_figure(fig, f"rwgttoS1_{name}", outdir=savedir)

    for name in tWLO_stpt1_hist_rwgt2: 
        h1 = tWLO_stpt1_hist_rwgt2[name]
        h2 = tWLO_stpt2_hist_rwgt2[name]

        fig, ax, rax = make_hist_ratio_plot_single(h1, h2, "Sample1", "Sample2", name, "Reweighted to Starting Point of Sample2")
        rax.set_ylabel('Ratio(S1/S2)')
        save_figure(fig, f"rwgttoS2_{name}", outdir=savedir)


    ##########################################
    # 1D Quad Plots: histEFT vs MGstandalone #
    ##########################################

    tWLO_stpt2_sow = "/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/TTto2L2Nu_stpt2_sow.pkl.gz"
    standAlone = plotTools.read_MGstandalone_txt("/users/hnelson2/afs/mc_production/standAloneMG/tW_SMEFTsim_top_xsec.txt")

    for wc in standAlone.keys():
        stpt2_list = get_rwgt_xsec_list(file=tWLO_stpt2_sow, wc_name=wc, wc_max=10.0)
        standalone_list = get_scatter_xsec_list(standAlone[wc], 5.626, 0.02046)

        samples = {'tWEFT':{'data':stpt2_list, 'type':'plot', 'label':'EFT sample'}, 
                    'MG':{'data':standalone_list, 'type':'scatter', 'label':'MG'}}
        fig, ax = make_1D_quad_plot(samples, wc)

        save_figure(fig, f"1dquad_{wc}", outdir=savedir)


    ##################################$
    # 1D Quad Plots: with LHEmmnl cut #
    ##################################$

    tWLO_stpt2_LHEmmnl_sow = "/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_LHEmmnl_sow.pkl.gz"
    standAlone_LHEmmnl = plotTools.read_MGstandalone_txt("/users/hnelson2/afs/mc_production/standAloneMG/tW_LHEcut_xsec.txt")

    for wc in standAlone_LHEmmnl.keys():
        stpt2_list = get_rwgt_xsec_list(file=tWLO_stpt2_LHEmmnl_sow, wc_name=wc, wc_max=12.0, hist_name="sow_mmnl_100")    
        #the factor of .105 is included because the MG standalone xsec was calculated without forcing the leptonic decay
        # 10.5% is the branching fraction to dilep final state
        standalone_list = get_scatter_xsec_list(standAlone_LHEmmnl[wc], 0.2877*0.105, 0.001149*0.105)  

        samples = {'tWEFT':{'data':stpt2_list, 'type':'plot', 'label':'EFT sample'}, 
                    'MG':{'data':standalone_list, 'type':'scatter', 'label':'MG'}}

        fig, ax = make_1D_quad_plot(samples, wc, const=0.105)
        ax.set_title("LHE.mmnl>100")
        ax.set_xlim([-8, 8])
        if 'cleQt1' in wc: 
            ax.set_ylim([0, 2.0])
        elif 'cleQt3' in wc: 
            ax.set_ylim([0, 8.0])
        elif 'cQl3' in wc: 
            ax.set_ylim([0, 5])

        save_figure(fig, f"1dquad_LHEmmnl100_{wc}", outdir=savedir)


    #########################################
    # 1D Quad Plots: dim6top vs SMEFTsimtop #
    #########################################

    SMEFTsim_standAlone = plotTools.read_MGstandalone_txt("/users/hnelson2/afs/mc_production/standAloneMG/tW_LHEcut_xsec.txt")
    dim6top_standAlone = plotTools.read_MGstandalone_txt("/users/hnelson2/afs/mc_production/standAloneMG/tW_dim6top_xsec.txt")

    wc_conversion = {'cQl311':'cQl31', 'cQl322':'cQl32', 'cQl333':'cQl33', 
                 'cleQt3Re11':'ctlT1', 'cleQt3Re22':'ctlT2', 'cleQt3Re33':'ctlT3', 
                 'cleQt1Re11':'ctlS1', 'cleQt1Re22':'ctlS2', 'cleQt1Re33':'ctlS3'}

    for wc in wc_conversion: 
        dim6wc = wc_conversion[wc]
        SMEFTsim_lst = get_scatter_xsec_list(SMEFTsim_standAlone[wc], 0.2877*0.105, 0.001149*0.105)
        dim6top_lst = get_scatter_xsec_list(dim6top_standAlone[dim6wc], 0.2483*0.105, 0.0006718*0.105)   

        samples = {'SMEFTsim':{'data':SMEFTsim_lst, 'type':'scatter', 'label':'SMEFTsim top'}, 
                    'dim6top':{'data':dim6top_lst, 'type':'scatter', 'label':'dim6top'}}

        fig, ax = make_1D_quad_plot(samples, wc, const=0.105)
        ax.set_title("LHE.mmnl>100")
        ax.set_xlim([-8, 8])
        if 'cleQt1' in wc: 
            ax.set_ylim([0, 1.6])
        elif 'cleQt3' in wc: 
            ax.set_ylim([0, 15])
        elif 'cQl3' in wc: 
            ax.set_ylim([0, 8])

        save_figure(fig, f"1dquad_dim6vsSMEFTsim_{wc}", outdir=savedir)


    #########################################
    # 2 EFT Samples Reweighted to eachother # 
    #########################################

    tWLO_stpt1_histEFT = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt1_histEFT.pkl.gz", allow_empty=False)
    tWLO_stpt2_histEFT = utils.get_hist_from_pkl("/users/hnelson2/afs/ttbarEFT/analysis/mc_validation/tWLO_stpt2_histEFT.pkl.gz", allow_empty=False)
 
    # starting point of tWLO_rwgt1 sample
    rwgt1 = {"ctGIm": -0.2, "ctGRe":-0.2, "cHQ3": 1.5, "ctWRe": -1.0, "cbWRe": -10.0, "cHtbRe": 4.0,
         "cQl311": 4.0, "cQl322": 4.1, "cQl333": 4.2,
         "cleQt3Re11": 12.0, "cleQt3Re22": 12.1, "cleQt3Re33": 12.2,
         "cleQt1Re11": 15.0, "cleQt1Re22": 15.1, "cleQt1Re33": 15.2}

    # starting point of tWLO_rwgt2 sample
    rwgt2 = {"ctGIm": -0.5, "ctGRe":-0.5, "cHQ3": 1.5, "ctWRe": -1.0, "cbWRe": -8.0, "cHtbRe": 4.0,
         "cQl311": 8.0, "cQl322": 8.0, "cQl333": 8.0,
         "cleQt3Re11": 12.0, "cleQt3Re22": 12.0, "cleQt3Re33": 12.0,
         "cleQt1Re11": 15.0, "cleQt1Re22": 15.0, "cleQt1Re33": 15.0} 

    # not including ratio
    for name in tWLO_stpt1_histEFT:
        samples={'stpt1': {'h':tWLO_stpt1_histEFT[name], 'rwgt': rwgt1, 
                    'plot_options': {'stack':False, 'yerr':False, 'linewidth':2, 'label':'sample1'}},
                'stpt2': {'h':tWLO_stpt2_histEFT[name], 'rwgt': rwgt1, 
                    'plot_options': {'stack':False, 'yerr':False, 'linewidth':2, 'label':'sample2'}}}

        ax_options = {'xlabel':name, 'ylabel':'Events', 'title': 'Reweighted to starting point of stpt1'}
        fig, ax = make_histEFT_plot_multi(samples, **ax_options)
        save_figure(fig, f"EFTrwgt1_{name}", outdir=savedir)

        samples={'stpt1': {'h':tWLO_stpt1_histEFT[name], 'rwgt': rwgt2, 
                    'plot_options': {'stack':False, 'yerr':False, 'linewidth':2, 'label':'sample1'}},
                'stpt2': {'h':tWLO_stpt2_histEFT[name], 'rwgt': rwgt2, 
                    'plot_options': {'stack':False, 'yerr':False, 'linewidth':2, 'label':'sample2'}}}
        ax_options = {'xlabel':name, 'ylabel':'Events', 'title': 'Reweighted to starting point of stpt2'}
        fig, ax = make_histEFT_plot_multi(samples, **ax_options)
        save_figure(fig, f"EFTrwgt2_{name}", outdir=savedir)

    # including ratio
    for name in tWLO_stpt1_histEFT:
        h_denom={'h':tWLO_stpt1_histEFT[name].as_hist(rwgt1), 
                    'plot_options': {'histtype': 'step', 'stack':False, 'yerr':False, 'linewidth':2, 'label':'stpt1'}} 
        samples={'stpt2': {'h':tWLO_stpt2_histEFT[name].as_hist(rwgt1),
                    'plot_options': {'histtype': 'step', 'stack':False, 'yerr':False, 'linewidth':2, 'label':'stpt2'}}}
        
        ax_options = {'xlabel':name, 'ylabel':'Events', 'title': 'Reweighted to starting point of stpt1'}
        fig, ax, rax = make_histEFT_ratio_plot_multi(samples, h_denom)
        rax.set_ylim([0,2])
        ax.set(**ax_options)
        save_figure(fig, f"EFTrwgt1_ratio_{name}", outdir=savedir)

        h_denom={'h':tWLO_stpt2_histEFT[name].as_hist(rwgt2), 
                    'plot_options': {'histtype': 'step', 'stack':False, 'yerr':False, 'linewidth':2, 'label':'stpt1'}} 
        samples={'stpt2': {'h':tWLO_stpt2_histEFT[name].as_hist(rwgt2), 
                    'plot_options': {'histtype': 'step', 'stack':False, 'yerr':False, 'linewidth':2, 'label':'stpt2'}}}
        
        ax_options = {'xlabel':name, 'ylabel':'Events', 'title': 'Reweighted to starting point of stpt2'}
        fig, ax, rax = make_histEFT_ratio_plot_multi(samples, h_denom)
        rax.set_ylim([0,2])
        ax.set(**ax_options)
        save_figure(fig, f"EFTrwgt2_ratio_{name}", outdir=savedir)

if __name__ == "__main__":
    main()