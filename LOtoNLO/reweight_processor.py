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
from hist import Hist
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
        
        # pt_bins = [*list(range(0, 500, 10)), *list(range(500, 600, 20)), *list(range(600, 800, 50)), *list(range(800, 1000, 100)), *list(range(1000, 2001, 500))]

        # much better for the ratio plot's smoothness
        pt_bins = [*list(range(0, 300, 10)), *list(range(300, 400, 20)), *list(range(400, 800, 50)), *list(range(800, 1000, 100)), *list(range(1000, 2001, 500))]


        pt_bins = [*list(range(0, 300, 10)), *list(range(300, 400, 20)), *list(range(400, 650, 50)), 650, 750, 850, 1000, 1500] #*list(range(700, 1000, 100)), *list(range(1000, 2001, 500))]

        self._histo_dict = {
            # "avg_toppt":  Hist(hist.axis.Regular(bins=200, start=0, stop=2000, name="avg_toppt", label="average $p_T$ of the top quarks [GeV]"), storage="weight"),
            "avg_toppt":  Hist(hist.axis.Variable(pt_bins, name="avg_toppt", label="average $p_T$ of the top quarks [GeV]"), storage="weight"),
            "nEvents":    Hist(hist.axis.Regular(bins=2, start=0, stop=1, name='nEvents', label="number of events"))
        }


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


        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None

        genpart = events.GenPart
        is_final_mask = genpart.hasFlags(["fromHardProcess","isLastCopy"])

        ######## Top selection ########

        gen_top = ak.pad_none(genpart[is_final_mask & (abs(genpart.pdgId) == 6)],2)
        avg_toppt = np.divide((gen_top[:,0].pt + gen_top[:,1].pt), 2)

        ######## Normalization ########
        if eft_coeffs is None:
            genw = events["genWeight"]
        else:
            genw = np.ones_like(events['event'])

        # genw = events['genWeight']
        counts = np.ones_like(events['event'])
        # else:                                       # If this is not an eft sample, get the genWeight
        #     genw = np.ones_like(events['event'])

        norm = genw*(xsec/sow)


        ######## Fill histos ########
        hout = self._histo_dict

        avg_toppt_info = {
            "avg_toppt" : avg_toppt,
            "weight"    : norm,
        }

        nEvents_info = {
            "nEvents" : counts,
        }

        hout['avg_toppt'].fill(**avg_toppt_info)
        hout['nEvents'].fill(**nEvents_info)

        # for var_name, var_values in variables_to_fill.items():
        #     fill_info = {
        #         var_name    : var_values,
        #         "weight"    : norm,
        #     }

        #     hout[var_name].fill(**fill_info)

        return hout


    def postprocess(self, accumulator):
        return accumulator
