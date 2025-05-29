#!/usr/bin/env python
import numpy as np
import awkward as ak

np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor
from coffea.analysis_tools import PackedSelection

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = False

import hist
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth

# Get the lumi for the given year
def get_lumi(year):
    lumi_dict = {
        "2016APV": 19.52,
        "2016": 16.81,
        "2017": 41.48,
        "2018": 59.83
    }
    if year not in lumi_dict.keys():
        raise Exception(f"(ERROR: Unknown year \"{year}\".")
    else:
        return(lumi_dict[year])

# Clean the objects
def is_clean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)

# Main analysis processor
class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst = None, dtype=np.float32, do_errors=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst

        self._dtype = dtype
        self._do_errors = do_errors

        print("\n\n")
        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)
        print("\n\n")

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)

        self._histo_dict = {
            "sow" :     HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow", label="sum of weights for all events"), 
                            wc_names=wc_names_lst, 
                            label="Events"),
            "sow_norm": HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow_norm", label="normalized sum of weights for all events"), 
                            wc_names=wc_names_lst, 
                            label="Events"), 
            "nEvents":  HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="nEvents", label="number of events"), 
                            wc_names=wc_names_lst, 
                            label="Events"),
            "lhe_mmnl":   HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=20, start=0, stop=200, name='lhe_mmnl', label='LHE invariant mass of initial leptons[GeV]'),
                            wc_names=wc_names_lst,
                            label="Events"),
            "lhe_mmnl_100":   HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=20, start=0, stop=200, name='lhe_mmnl_100', label='LHE invariant mass of initial leptons[GeV]'),
                            wc_names=wc_names_lst,
                            label="Events"),
            "sow_mmnl_100" :HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow_mmnl_100", label="sum of weights for events mmnl>100"), 
                            wc_names=wc_names_lst, 
                            label="Events"),
        }

        # Set the list of hists to to fill
        if hist_lst is None:
            self._hist_lst = list(self._histo_dict.keys())
        else:
            for h in hist_lst:
                if h not in self._histo_dict.keys():
                    raise Exception(f"Error: Cannot specify hist \"{h}\", it is not defined in self._histo_dict")
            self._hist_lst = hist_lst

        print("hist_lst: ", self._hist_lst)


    @property
    def columns(self):
        return self._columns

    def process(self, events):

        # Dataset parameters
        dataset = events.metadata['dataset']

        isData = self._samples[dataset]["isData"]
        hist_axis_name = self._samples[dataset]["histAxisName"]
        year   = self._samples[dataset]['year']
        xsec   = self._samples[dataset]['xsec']
        sow    = self._samples[dataset]['nSumOfWeights']

        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None
        # eft_w2_coeffs = efth.calc_w2_coeffs(eft_coeffs,self._dtype) if (self._do_errors and eft_coeffs is not None) else None


        # Initialize objects
        lhepart = events.LHEPart

        genpart = events.GenPart

        ######## mmnl cut #######
        '''invariant mass of the initial leptons (not from the top decay)
        the tW samples were produced using the MG process definitions: 
          p p > t l- vl~, (t > l+ vl b) 
          p p > t~ l+ vl, (t~ > l- vl~ b~
        when this was used to make the nanogen samples, lhepart[:,5]+lhepart[:,6] are always the lv pair NOT from the top decay 
        confirmation of this can be found in ttbarEFT/analysis/mc_validation/check_tW_LHEordering.ipynb
        '''
        lhe_mmnl = (lhepart[:,5]+lhepart[:,6]).mass

        # mask of just inv mass > 100 to match 100 = mmnl used in MG standalone xsec calculation
        mmnl_more100 = ak.fill_none(lhe_mmnl>100, False)

        selections = PackedSelection()
        selections.add('mmnl', mmnl_more100)
        event_selection_mask = selections.all('mmnl')

        ######## Normalization ########

        norm = (xsec/sow)

        if eft_coeffs is None:
            genw = events["genWeight"]
        else:
            genw = np.ones_like(events['event'])

        counts = np.ones_like(events['event'])
        event_weights = genw*norm

        ######## Fill histos ########

        hout = self._histo_dict

        sow_fill_info = {
            "sow"       : counts, 
            "process"   : hist_axis_name, 
            "weight"    : genw, 
            "eft_coeff" : eft_coeffs,
        }

        sow_norm_fill_info = {
            "sow_norm"  : counts, 
            "process"   : hist_axis_name, 
            "weight"    : event_weights*norm, 
            "eft_coeff" : eft_coeffs,
        }

        sow_mmnl_100_fill_info = {
            "sow_mmnl_100"  : counts[event_selection_mask],
            "process"       : hist_axis_name,
            "weight"        : event_weights[event_selection_mask],
            "eft_coeff"     : eft_coeffs[event_selection_mask],
        }

        lhe_mmnl_fill_info = {
            "lhe_mmnl" : lhe_mmnl,
            "process"   : hist_axis_name,
            "weight"    : event_weights,
            "eft_coeff" : eft_coeffs,
        }

        lhe_mmnl100_fill_info = {
            "lhe_mmnl_100" : lhe_mmnl[event_selection_mask],
            "process"   : hist_axis_name,
            "weight"    : event_weights[event_selection_mask],
            "eft_coeff" : eft_coeffs[event_selection_mask],
        }

        # Here, weight = counts instead of wgts because this hist is just counting the number
        # of raw events in the file, not effected by the weighting of the event
        nevents_fill_info = {
            "nEvents"   : counts[event_selection_mask], 
            "process"   : hist_axis_name, 
            "weight"    : counts[event_selection_mask],
            "eft_coeff" : None,
        }

        hout['nEvents'].fill(**nevents_fill_info)
        hout['sow'].fill(**sow_fill_info)
        hout['sow_norm'].fill(**sow_norm_fill_info)
        hout['sow_mmnl_100'].fill(**sow_mmnl_100_fill_info)
        hout['lhe_mmnl'].fill(**lhe_mmnl_fill_info)
        hout['lhe_mmnl_100'].fill(**lhe_mmnl100_fill_info)

        return hout


    def postprocess(self, accumulator):
        return accumulator
