import numpy as np
import awkward as ak
np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = False

import hist
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth
import ttbarEFT.modules.corrections as tt_cor 
import topcoffea.modules.corrections as tc_cor
from ttbarEFT.modules.processor_tools import calc_eft_weights

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

        # Create the histograms with new scikit hist
        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        self._histo_dict = {
            "djr_10_all": HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=80, start=0, stop=4, name="djr_10_all", label="DJR  0 to 1"), 
                            wc_names=wc_names_lst, 
                            label="Events"), 
            "djr_10_0p": HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=80, start=0, stop=4, name="djr_10_0p", label="DJR  0 to 1"), 
                            wc_names=wc_names_lst, 
                            label="Events"), 
            "djr_10_1p": HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=80, start=0, stop=4, name="djr_10_1p", label="DJR  0 to 1"), 
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
        hist_axis_name = self._samples[dataset]["histAxisName"]

        year   = self._samples[dataset]['year']
        xsec   = self._samples[dataset]['xsec']
        sow    = self._samples[dataset]['nSumOfWeights']

        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None
        eft_w2_coeffs = efth.calc_w2_coeffs(eft_coeffs,self._dtype) if (self._do_errors and eft_coeffs is not None) else None

        djr_10 = events.LHEWeight.DJR10
        nMEPartons = events.Generator.nMEPartonsFiltered

        djr_10_0p = djr_10[nMEPartons == 0]
        djr_10_1p = djr_10[nMEPartons == 1]

        wgts = np.ones_like(events['event'])
        if eft_coeffs is None:
            # Basically any central MC samples
            wgts = events["genWeight"]

        ####### Fill Histogram #######
        hout = self._histo_dict

        djr_10_fill_info = {
            "djr_10_all"  	: np.log10(djr_10),
            "process"    	: hist_axis_name,
            "weight"    	: wgts,
            "eft_coeff" 	: eft_coeffs,
        }

        djr_10_0p_fill_info = {
            "djr_10_0p"  	: np.log10(djr_10_0p),
            "process"    	: hist_axis_name,
            "weight"    	: wgts[nMEPartons == 0],
            "eft_coeff" 	: eft_coeffs[nMEPartons == 0],
        }

        djr_10_1p_fill_info = {
            "djr_10_1p"  	: np.log10(djr_10_1p),
            "process"    	: hist_axis_name,
            "weight"    	: wgts[nMEPartons == 1],
            "eft_coeff" 	: eft_coeffs[nMEPartons == 1],
        }

        hout["djr_10_all"].fill(**djr_10_fill_info)
        hout["djr_10_0p"].fill(**djr_10_0p_fill_info)
        hout["djr_10_1p"].fill(**djr_10_1p_fill_info)

        return hout


    def postprocess(self, accumulator):
        return accumulator





