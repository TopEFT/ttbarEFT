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
            "nEvents":  HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="nEvents", label="number of events"), 
                            wc_names=wc_names_lst, 
                            label="Events"),
            "gen_mtt":   HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=65, start=0, stop=1300, name='gen_mtt', label='GEN invariant mass of tops [GeV]'),
                            wc_names=wc_names_lst,
                            label="Events"),
            "lhe_mtt":   HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=65, start=0, stop=1300, name='lhe_mtt', label='LHE invariant mass of tops [GeV]'),
                            wc_names=wc_names_lst,
                            label="Events"),
            "sow_0_700" :HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow_0_700", label="sum of weights for all events"), 
                            wc_names=wc_names_lst, 
                            label="Events"),
            "sow_700_900" :HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow_700_900", label="sum of weights for all events"), 
                            wc_names=wc_names_lst, 
                            label="Events"),
            "sow_900_Inf" :HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow_900_Inf", label="sum of weights for all events"), 
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
        is_final_mask = genpart.hasFlags(["fromHardProcess","isLastCopy"])

        ######## Top selection ########

        # mttbar calculation using gen particles 
        gen_top = ak.pad_none(genpart[is_final_mask & (abs(genpart.pdgId) == 6)],2)
        gen_mtt = (gen_top[:,0] + gen_top[:,1]).mass

        # mttbar calculation using LHE particles. This is identical to the filter used in the central sample production
        lhe_mtt = (lhepart[:,0]+lhepart[:,1]+lhepart[:,2]+lhepart[:,3]+lhepart[:,4]+lhepart[:,5]+lhepart[:,6]+lhepart[:,7]).mass

        mtt_less700 = ak.fill_none(lhe_mtt<700, False)
        mtt_more700 = ak.fill_none(lhe_mtt >=700, False)
        mtt_less900 = ak.fill_none(lhe_mtt <=900, False)
        mtt_more900 = ak.fill_none(lhe_mtt > 900, False)

        mtt_0_700 = PackedSelection()
        mtt_0_700.add('less700', mtt_less700)

        mtt_700_900 = PackedSelection()
        mtt_700_900.add('more700', mtt_more700)
        mtt_700_900.add('less900', mtt_less900)

        mtt_900_Inf = PackedSelection()
        mtt_900_Inf.add('more900', mtt_more900)

        mask_0_700 = mtt_0_700.all('less700')
        mask_700_900 = mtt_700_900.all('more700', 'less900')
        mask_900_Inf = mtt_900_Inf.all('more900')

        ######## Normalization ########

        norm = 1

        if eft_coeffs is None:
            genw = events["genWeight"]
        else:
            genw = np.ones_like(events['event'])

        event_weights = norm*genw

        counts = np.ones_like(events['event'])

        ######## Fill histos ########

        hout = self._histo_dict

        variables_to_fill = {
            "sow"     : counts,
            "gen_mtt"  : gen_mtt,
            "lhe_mtt" : lhe_mtt,
        }

        for var_name, var_values in variables_to_fill.items():
            if var_name not in self._hist_lst:
                print(f"Skipping \"{var_name}\", it is not in the list of hists to include")
                continue

            fill_info = {
                var_name    : var_values,
                "process"   : hist_axis_name,
                "weight"    : event_weights,
                "eft_coeff" : eft_coeffs,
            }

            hout[var_name].fill(**fill_info)

        # fill_nevents = {
        #     "nEvents"   : counts,
        #     "process"   : hist_axis_name,
        #     "weight"    : event_weights,
        #     "eft_coeff" : None,
        # }

        # fill_sow_0_700 = {
        #     "sow_0_700" : counts[mask_0_700],
        #     "process"   : hist_axis_name,
        #     "weight"    : event_weights[mask_0_700],
        #     "eft_coeff" : eft_coeffs[mask_0_700],
        # }

        # fill_sow_700_900 = {
        #     "sow_700_900" : counts[mask_700_900],
        #     "process"   : hist_axis_name,
        #     "weight"    : event_weights[mask_700_900],
        #     "eft_coeff" : eft_coeffs[mask_700_900],
        # }

        # fill_sow_900_Inf = {
        #     "sow_900_Inf" : counts[mask_900_Inf],
        #     "process"   : hist_axis_name,
        #     "weight"    : event_weights[mask_900_Inf],
        #     "eft_coeff" : eft_coeffs[mask_900_Inf],
        # }

        # hout['nEvents'].fill(**fill_nevents)
        # hout['sow_0_700'].fill(**fill_sow_0_700)
        # hout['sow_700_900'].fill(**fill_sow_700_900)
        # hout['sow_900_Inf'].fill(**fill_sow_900_Inf)

        return hout


    def postprocess(self, accumulator):
        return accumulator
