import os 
import gzip
import cloudpickle

import awkward as ak
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import mplhep as hep

import hist
from hist import Hist

# from coffea.analysis_tools import PackedSelection
from topcoffea.modules import utils
import topcoffea.modules.eft_helper as efth

plot_axes = {
    'ljets': 'number of light jets',
    'mtt': 'mtt [GeV]',
    'pttt': r'$p_T$(tt) [GeV]',
    'top1pt': r'leading top $p_T$ [GeV]',
    'top2pt': r'subleading top $p_T$ [GeV]',
    'top1eta': r'leading top $\eta$',
    'top2eta': r'subleading top $\eta$',
    'top1phi': r'leading top $\phi$',
    'top2phi': r'subleading top $\phi$',
    'mll': 'mll [GeV]',
    'ptll': r'$p_T$(ll) [GeV]',
    'l0pt': r'leading lepton $p_T$ [GeV]',
    'l1pt': r'subleading lepton $p_T$ [GeV]',
    'l0eta': r'leading lepton $\eta$',
    'l1eta': r'subleading lepton $\eta$',
    'l0phi': r'leading lepton $\phi$',
    'l1phi': r'subleading lepton $\phi$',
    'dr_leps': r'$\Delta$R(leading, subleading leptons)',
    'j0pt': r'leading jet $p_T$',
    'j0eta': r'leading jet $\eta$',
    'j0phi': r'leadin gjet $\phi$',
}

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

    yvals_num = num_hist.values()
    yvals_denom = denom_hist.values()
    sigma_num = np.sqrt(num_hist.variances())
    sigma_denom = np.sqrt(denom_hist.variances())

    ratio = np.divide(yvals_num, yvals_denom)

    # calculation for error propagation for ratio = yavls_num/yvals_denom
    # generally, z=x/y; sigma_z = abs(z)sqrt((sigma_x/x)^2+(sigma_y/y)^2)
    sigma_y = np.multiply(np.abs(ratio), np.sqrt(np.add(np.square(np.divide(sigma_num, yvals_num)), np.square(np.divide(sigma_denom, yvals_denom)))))

    return sigma_y


def save_figure(fig, figname, outdir=''):
    outname=os.path.join(outdir, f"{figname}.png")
    fig.savefig(outname, bbox_inches='tight')
    print(f'plot saved to {outname}')
    plt.close(fig)


def plot_single(hist):
    hep.style.use('CMS')
    fig, ax = plt.subplots()

    hep.histplot(hist, ax=ax, stack=False, yerr=True, linewidth=2, label='new SMEFTsim')
    ax.legend()

    return fig, ax 

def multi_EFT_plot(h_SM, h_EFT1, h_EFT2):
    hep.style.use('CMS')
    fig, ax = plt.subplots()

    hep.histplot(h_SM, ax=ax, stack=False, yerr=True, linewidth=3, label='no EFT', color='black')
    hep.histplot(h_EFT1, ax=ax, stack=False, yerr=True, linewidth=3, label='EFT option 1')
    hep.histplot(h_EFT2, ax=ax, stack=False, yerr=True, linewidth=3, label='EFT option 2')

    ax.legend(loc='upper right')

    return fig, ax


def plot_comparison(h_smeft_old, h_smeft_new, h_powheg):

    bin_widths = h_powheg.axes[0].edges
    centers = h_powheg.axes[0].centers
    val_denom = h_powheg.values()

    ratio_smeft_old = (h_smeft_old.values())/val_denom
    ratio_smeft_new = (h_smeft_new.values())/val_denom

    hep.style.use('CMS')
    # Initialize figure and axes
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)

    hep.histplot(h_powheg, ax=ax, stack=False, yerr=True, linewidth=2, label='powheg', color='black')
    hep.histplot(h_smeft_old, ax=ax, stack=False, yerr=True, linewidth=2, label='old EFT')
    hep.histplot(h_smeft_new, ax=ax, stack=False, yerr=True, linewidth=2, label='new EFT')

    rax.scatter(centers, ratio_smeft_old, label='old EFT')
    rax.scatter(centers, ratio_smeft_new, label='new EFT')

    ax.legend(loc='upper right')
    ax.ticklabel_format(style="sci", scilimits=(-3, 3), useMathText=True)
    rax.axhline(y=1.0, color='gray', linestyle='--')

    return fig, ax, rax

def plot_powheg_MLM_comparison(h_smeft, h_powheg, h_MLM, SMEFTsim_label='SMEFTsim', SMEFTsim_color='red'):

    bin_widths = h_powheg.axes[0].edges
    centers = h_powheg.axes[0].centers
    val_denom = h_powheg.values()

    ratio_mg_MLM = (h_MLM.values())/val_denom
    ratio_smeft = (h_smeft.values())/val_denom

    uncert_ratio_mgMLM = get_ratio_uncertainty(h_MLM, h_powheg)
    uncert_ratio_smeft = get_ratio_uncertainty(h_smeft, h_powheg)

    hep.style.use('CMS')
    # Initialize figure and axes
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)

    nom_MG_color = '#63CCA3'

    hep.histplot(h_powheg, ax=ax, stack=False, yerr=True, linewidth=2, label='Powheg NLO', color='black')
    hep.histplot(h_MLM, ax=ax, stack=False, yerr=True, linewidth=2, label='MadGraph MLM', color=nom_MG_color)
    hep.histplot(h_smeft, ax=ax, stack=False, yerr=True, linewidth=2, label=SMEFTsim_label, color=SMEFTsim_color)

    # rax.scatter(centers, ratio_mg_MLM, label='MadGraph MLM')
    # rax.scatter(centers, ratio_smeft, label='SMEFTsim')

    rax.scatter(centers, ratio_mg_MLM, color=nom_MG_color, label='MadGraph MLM')
    rax.errorbar(centers, ratio_mg_MLM, xerr = None, yerr = uncert_ratio_mgMLM, capsize=5, ls='none', color=nom_MG_color)

    rax.scatter(centers, ratio_smeft, color=SMEFTsim_color, label=SMEFTsim_label)
    rax.errorbar(centers, ratio_smeft, xerr=None, yerr=uncert_ratio_smeft, capsize=5, ls='none', color=SMEFTsim_color)

    ax.legend(loc='best')
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)
    rax.axhline(y=1.0, color='gray', linestyle='--')

    return fig, ax, rax

# def plot_powheg_MLMcomparison(h_smeft_old, h_smeft_new, h_powheg, h_MLM):

#     bin_widths = h_powheg.axes[0].edges
#     centers = h_powheg.axes[0].centers
#     val_denom = h_powheg.values()

#     ratio_mg_MLM = (h_MLM.values())/val_denom
#     ratio_smeft_old = (h_smeft_old.values())/val_denom
#     ratio_smeft_new = (h_smeft_new.values())/val_denom

#     hep.style.use('CMS')
#     # Initialize figure and axes
#     fig, (ax, rax) = plt.subplots(
#         nrows=2,
#         ncols=1,
#         figsize=(10,12),
#         gridspec_kw={'height_ratios': (3, 1)},
#         sharex=True
#     )
#     fig.subplots_adjust(hspace=.1)

#     hep.histplot(h_powheg, ax=ax, stack=False, yerr=True, linewidth=2, label='Powheg', color='black')
#     hep.histplot(h_MLM, ax=ax, stack=False, yerr=True, linewidth=2, label='MadGraph MLM')
#     hep.histplot(h_smeft_old, ax=ax, stack=False, yerr=True, linewidth=2, label='old SMEFTsim')
#     hep.histplot(h_smeft_new, ax=ax, stack=False, yerr=True, linewidth=2, label='new SMEFTsim')

#     rax.scatter(centers, ratio_mg_MLM, label='MadGraph MLM')
#     rax.scatter(centers, ratio_smeft_old, label='old SMEFTsim')
#     rax.scatter(centers, ratio_smeft_new, label='new SMEFTsim')
    
#     ax.legend(loc='best')
#     ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)
#     rax.axhline(y=1.0, color='gray', linestyle='--')

#     return fig, ax, rax

def plot_EFT_comparison(h_smeft_new, h_SM_new, h_smeft_old, h_SM_old):

    bin_widths = h_smeft_new.axes[0].edges
    centers = h_smeft_new.axes[0].centers

    ratio_smeft_old = (h_smeft_old.values())/(h_SM_old.values())
    ratio_smeft_new = (h_smeft_new.values())/(h_SM_new.values())

    hep.style.use('CMS')
    # Initialize figure and axes
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)

    hep.histplot(h_SM_old, ax=ax, stack=False, yerr=True, linewidth=2, histtype='fill', color='#5790fc', label='old, SM')
    hep.histplot(h_SM_new, ax=ax, stack=False, yerr=True, linewidth=2, histtype='fill', color='#f89c20', label='new, SM')
    hep.histplot(h_smeft_old, ax=ax, stack=False, yerr=True, linewidth=2, histtype='step', color='blue', label='old, EFT rwgt')
    hep.histplot(h_smeft_new, ax=ax, stack=False, yerr=True, linewidth=2, histtype='step', color='#eb5b02', label='new, EFTrwgt')

    rax.scatter(centers, ratio_smeft_old, label='old EFT', color='blue')
    rax.scatter(centers, ratio_smeft_new, label='new EFT', color='#eb5b02')

    ax.legend(loc='best')
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)
    rax.set_ylabel('EFT/SM')
    rax.axhline(y=1.0, color='gray', linestyle='--')

    return fig, ax, rax


def plot_powheg_old_new_comparison(h_smeft_old, h_smeft_new, h_powheg):

    bin_widths = h_powheg.axes[0].edges
    centers = h_powheg.axes[0].centers
    val_denom = h_powheg.values()

    ratio_smeft_old = (h_smeft_old.values())/val_denom
    ratio_smeft_new = (h_smeft_new.values())/val_denom

    uncert_smeft_old = get_ratio_uncertainty(h_smeft_old, h_powheg)
    uncert_smeft_new = get_ratio_uncertainty(h_smeft_new, h_powheg)

    hep.style.use('CMS')
    # Initialize figure and axes
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)

    # nom_MG_color = '#63CCA3'
    new_smeft_color = 'magenta'
    old_smeft_color = 'blue'

    hep.histplot(h_powheg, ax=ax, stack=False, yerr=True, linewidth=2, label='Powheg NLO', color='black')
    hep.histplot(h_smeft_new, ax=ax, stack=False, yerr=True, linewidth=2, label='new SMEFTsim', color=new_smeft_color)
    hep.histplot(h_smeft_old, ax=ax, stack=False, yerr=True, linewidth=2, label='old SMEFTsim', color=old_smeft_color)

    # rax.scatter(centers, ratio_mg_MLM, label='MadGraph MLM')
    # rax.scatter(centers, ratio_smeft, label='SMEFTsim')

    rax.scatter(centers, ratio_smeft_new, color=new_smeft_color, label='new SMEFTsim')
    rax.errorbar(centers, ratio_smeft_new, xerr = None, yerr = uncert_smeft_new, capsize=5, ls='none', color=new_smeft_color)

    rax.scatter(centers, ratio_smeft_old, color=old_smeft_color, label='old SMEFTsim')
    rax.errorbar(centers, ratio_smeft_old, xerr=None, yerr=uncert_smeft_old, capsize=5, ls='none', color=old_smeft_color)

    ax.legend(loc='best')
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)
    rax.axhline(y=1.0, color='gray', linestyle='--')

    return fig, ax, rax


def plot_powheg_new_MLM_comparison(h_MLM, h_smeft_new, h_powheg):

    bin_widths = h_powheg.axes[0].edges
    centers = h_powheg.axes[0].centers
    val_denom = h_powheg.values()

    ratio_mg_MLM = (h_MLM.values())/val_denom
    ratio_smeft_new = (h_smeft_new.values())/val_denom

    uncert_ratio_mgMLM = get_ratio_uncertainty(h_MLM, h_powheg)
    uncert_ratio_smeft = get_ratio_uncertainty(h_smeft_new, h_powheg)

    hep.style.use('CMS')
    # Initialize figure and axes
    fig, (ax, rax) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(10,12),
        gridspec_kw={'height_ratios': (3, 1)},
        sharex=True
    )
    fig.subplots_adjust(hspace=.1)

    # nom_MG_color = '#63CCA3'
    new_smeft_color = 'magenta'
    MLM_color = '#63CCA3'

    hep.histplot(h_powheg, ax=ax, stack=False, yerr=True, linewidth=2, label='Powheg NLO', color='black')
    hep.histplot(h_smeft_new, ax=ax, stack=False, yerr=True, linewidth=2, label='new SMEFTsim', color=new_smeft_color)
    hep.histplot(h_MLM, ax=ax, stack=False, yerr=True, linewidth=2, label='MadGraph MLM', color=MLM_color)

    # rax.scatter(centers, ratio_mg_MLM, label='MadGraph MLM')
    # rax.scatter(centers, ratio_smeft, label='SMEFTsim')

    rax.scatter(centers, ratio_smeft_new, color=new_smeft_color, label='new SMEFTsim')
    rax.errorbar(centers, ratio_smeft_new, xerr = None, yerr = uncert_ratio_smeft, capsize=5, ls='none', color=new_smeft_color)

    rax.scatter(centers, ratio_mg_MLM, color=MLM_color, label='MadGraph MLM')
    rax.errorbar(centers, ratio_mg_MLM, xerr=None, yerr=uncert_ratio_mgMLM, capsize=5, ls='none', color=MLM_color)

    ax.legend(loc='best')
    ax.ticklabel_format(axis='y', style="sci", scilimits=(-3, 3), useMathText=True)
    rax.axhline(y=1.0, color='gray', linestyle='--')

    return fig, ax, rax


if __name__ == '__main__':
    
    f_new_smeft_SM = "SMpt_TT01j2lBSM_histos.pkl.gz"
    f_old_smeft_SM = "SMpt_oldSMEFTsim_histos.pkl.gz"
    f_powheg = "SMpt_powheg_histos.pkl.gz"
    f_mgMLM = "SMpt_mgMLM_histos.pkl.gz"

    f_new_smeft_EFT1 = "EFTpt1_TT01j2lBSM_histos.pkl.gz"
    f_new_smeft_EFT2 = "EFTpt2_TT01j2lBSM_histos.pkl.gz"

    f_old_smeft_EFT1 = "EFTpt1_oldSMEFTsim_histos.pkl.gz"
    f_old_smeft_EFT2 = "EFTpt2_oldSMEFTsim_histos.pkl.gz"

    new_smeft_SM_hists = utils.get_hist_from_pkl(f_new_smeft_SM, allow_empty=False)
    old_smeft_SM_hists = utils.get_hist_from_pkl(f_old_smeft_SM, allow_empty=False)
    powheg_hists = utils.get_hist_from_pkl(f_powheg, allow_empty=False)
    mgMLM_hists = utils.get_hist_from_pkl(f_mgMLM, allow_empty=False)

    new_smeft_EFT1_hists = utils.get_hist_from_pkl(f_new_smeft_EFT1, allow_empty=False)
    new_smeft_EFT2_hists = utils.get_hist_from_pkl(f_new_smeft_EFT2, allow_empty=False)

    old_smeft_EFT1_hists = utils.get_hist_from_pkl(f_old_smeft_EFT1, allow_empty=False)
    old_smeft_EFT2_hists = utils.get_hist_from_pkl(f_old_smeft_EFT2, allow_empty=False)

    ### SM comparison ###
    # keys_to_plot = new_smeft_hists.keys()
    # keys_to_plot = ['j0eta']
    # for name in keys_to_plot:
    #     fig, ax, rax = plot_comparison(
    #                         h_smeft_old=old_smeft_hists[name], 
    #                         h_smeft_new=new_smeft_hists[name], 
    #                         h_powheg=powheg_hists[name])

    #     ax.set_title("Reweighted to SM")
    #     ax.set_ylabel('Events')
    #     ax.set_xlabel('')
    #     if name in plot_axes.keys():
    #         rax.set_xlabel(plot_axes[name])
    #     else:
    #         rax.set_xlabel(name)

    #     save_figure(fig, name, "2610_plots")

    # print("mll", mgMLM_hists['mll'].values())
    # print("ptll", mgMLM_hists['ptll'].variances())

    # print(f"{mgMLM_hists['ptll'].variances()==None}")

    # print(f"{dir(mgMLM_hists['ptll'])}")

    ## Comparison between powheg, nominal madgraph, and SMEFTsim ###
    keys_to_plot = mgMLM_hists.keys()
    # keys_to_plot = ['njets']
    for name in keys_to_plot:
        if name == 'ptll': continue

        fig, ax, rax = plot_powheg_old_new_comparison(
                            h_smeft_old=old_smeft_SM_hists[name], 
                            h_smeft_new=new_smeft_SM_hists[name], 
                            h_powheg=powheg_hists[name])

        # fig, ax, rax = plot_powheg_MLM_comparison(
        #                     h_smeft=old_smeft_SM_hists[name], 
        #                     h_powheg=powheg_hists[name], 
        #                     h_MLM=mgMLM_hists[name], 
        #                     SMEFTsim_label = "old SMEFTsim", 
        #                     SMEFTsim_color='blue')

        ax.set_title("Reweighted to SM")
        ax.set_ylabel('Events')
        rax.set_ylabel('MG / Powheg')
        ax.set_xlabel('')
        if name in plot_axes.keys():
            rax.set_xlabel(plot_axes[name])
        else:
            rax.set_xlabel(name)
        rax.set_ylim([0.5, 1.5])

        save_figure(fig, f"PowhegOldNew_{name}", "301025_plots")

    # keys_to_plot = ['njets']
    for name in keys_to_plot:
        if name == 'ptll': continue

        fig, ax, rax = plot_powheg_new_MLM_comparison(
                            h_MLM=mgMLM_hists[name], 
                            h_smeft_new=new_smeft_SM_hists[name], 
                            h_powheg=powheg_hists[name])

        # fig, ax, rax = plot_powheg_MLM_comparison(
        #                     h_smeft=new_smeft_SM_hists[name], 
        #                     h_powheg=powheg_hists[name], 
        #                     h_MLM=mgMLM_hists[name], 
        #                     SMEFTsim_label = "new SMEFTsim", 
        #                     SMEFTsim_color='magenta')

        ax.set_title("Reweighted to SM")
        ax.set_ylabel('Events')
        rax.set_ylabel('MG / Powheg')
        ax.set_xlabel('')
        if name in plot_axes.keys():
            rax.set_xlabel(plot_axes[name])
        else:
            rax.set_xlabel(name)

        rax.set_ylim([0.5, 1.5])

        save_figure(fig, f"PowhegNewMLM_{name}", "301025_plots")

    ### Comparison of all 4 @ SM ###
    # keys_to_plot = new_smeft_SM_hists.keys()
    # # keys_to_plot = ['j0eta']
    # for name in keys_to_plot:
    #     fig, ax, rax = plot_powheg_MLMcomparison(
    #                         h_smeft_new=new_smeft_SM_hists[name], 
    #                         h_smeft_old = old_smeft_SM_hists[name],
    #                         h_powheg=powheg_hists[name], 
    #                         h_MLM=mgMLM_hists[name])

    #     ax.set_title("Reweighted to SM")
    #     ax.set_ylabel('Events')
    #     rax.set_ylabel('MG / Powheg')
    #     ax.set_xlabel('')
    #     if name in plot_axes.keys():
    #         rax.set_xlabel(plot_axes[name])
    #     else:
    #         rax.set_xlabel(name)

    #     save_figure(fig, f"all_{name}", "281025_plots")


    # ### EFT rwgt 1 ###
    # keys_to_plot = new_smeft_EFT1_hists.keys()

    # keys_to_plot = ['mtt']
    # for name in keys_to_plot:
    #     fig, ax, rax = plot_EFT_comparison(
    #                         h_smeft_new=new_smeft_EFT1_hists[name], 
    #                         h_SM_new=new_smeft_SM_hists[name], 
    #                         h_smeft_old=old_smeft_EFT1_hists[name], 
    #                         h_SM_old=old_smeft_SM_hists[name])

    #     ax.set_title("Reweighted to SM and EFT1")
    #     ax.set_ylabel('Events')
    #     ax.set_xlabel('')
    #     if name in plot_axes.keys():
    #         rax.set_xlabel(plot_axes[name])
    #     else:
    #         rax.set_xlabel(name)

    #     save_figure(fig, f'{name}_EFT1', '2610_plots')


    # ### EFT rwgt 2 ###
    # keys_to_plot = new_smeft_EFT2_hists.keys()

    # keys_to_plot = ['mtt']
    # for name in keys_to_plot:
    #     fig, ax, rax = plot_EFT_comparison(
    #                         h_smeft_new=new_smeft_EFT2_hists[name], 
    #                         h_SM_new=new_smeft_SM_hists[name], 
    #                         h_smeft_old=old_smeft_EFT2_hists[name], 
    #                         h_SM_old=old_smeft_SM_hists[name])

    #     ax.set_title("Reweighted to SM and EFT2")
    #     ax.set_ylabel('Events')
    #     ax.set_xlabel('')
    #     if name in plot_axes.keys():
    #         rax.set_xlabel(plot_axes[name])
    #     else:
    #         rax.set_xlabel(name)

    #     save_figure(fig, f'{name}_EFT2', '2610_plots')


    ### EFT comparison plots ###
    # hists_noEFT = utils.get_hist_from_pkl("2410_old_SMEFTsim_histos.pkl.gz", allow_empty=False)
    # hists_EFT1 = utils.get_hist_from_pkl("EFTpt1_SMEFTsim_histos.pkl.gz", allow_empty=False)
    # hists_EFT2 = utils.get_hist_from_pkl("EFTpt2_SMEFTsim_histos.pkl.gz", allow_empty=False)

    # fig, ax = multi_EFT_plot(
    #     h_SM=hists_noEFT['top1pt'], 
    #     h_EFT1=hists_EFT1['top1pt'], 
    #     h_EFT2=hists_EFT2['top1pt'])

    # ax.set_title("Top Quark Momentum")
    # ax.set_ylabel("Events")

    # save_figure(fig, 'top1pt', "presentation_plots")





