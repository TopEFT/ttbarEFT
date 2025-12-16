#!/usr/bin/env python
import copy
import coffea
import numpy as np
import awkward as ak
import json
import hist
import yaml

# from mt2 import mt2

from coffea import processor
# from coffea.analysis_tools import PackedSelection
# from coffea.nanoevents.methods import vector
# from coffea.lumi_tools import LumiMask

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from topcoffea.modules.paths import topcoffea_path
from topcoffea.modules.histEFT import HistEFT
# import topcoffea.modules.eft_helper as efth
# import topcoffea.modules.corrections as tc_cor
# import topcoffea.modules.event_selection as tc_es

# from ttbarEFT.modules.paths import ttbarEFT_path
# from ttbarEFT.modules.analysis_tools import make_mt2, get_lumi
# import ttbarEFT.modules.object_selection as tt_os
# import ttbarEFT.modules.event_selection as tt_es
# from ttbarEFT.modules.corrections import AttachElectronSF, AttachMuonSF, AttachMuonTrigSF, ApplyMuonPtCorr, ApplyJetVetoMaps

# from topcoffea.modules.get_param_from_jsons import GetParam

# get_tc_param = GetParam(topcoffea_path("params/params.json"))
# get_tt_param = GetParam(ttbarEFT_path("params/params.json"))

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, lep_cat, wc_names_lst=[], hist_lst=None, dtype=np.float32):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 
        self._lep_cat = lep_cat

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)

        histograms = {}
        histograms['njets'] = HistEFT(
            proc_axis,
            hist.axis.Regular(8, 0, 8, name='njets', label='njets'),
            wc_names = wc_names_lst,
            label=r'Events',
        )

        self._accumulator = histograms

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    def process(self, events):

        # Dataset parameters
        dataset         = events.metadata['dataset']
        isData          = self._samples[dataset]['isData']
        histAxisName    = self._samples[dataset]['histAxisName']
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']


        isEFT = hasattr(events, 'EFTfitCoefficients')    
        assert not (isEFT and isData), f"isEFT and isData cannot both be True. Check input samples."

        lep_cat = self._lep_cat

        datasets = ["Muon", "SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG"]
        for d in datasets:
            if dataset.startswith(d):
                dataset = dataset.split('_')[0]


        ######### EFT coefficients ##########
        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, 'EFTfitCoefficients') else None
        if eft_coeffs is not None:
            # Check to see if the ordering of WCs for this sample matches what want
            if self._samples[dataset]['WCnames'] != self._wc_names_lst:
                eft_coeffs = efth.remap_coeffs(self._samples[dataset]['WCnames'], self._wc_names_lst, eft_coeffs)

        # Initialize the out object
        hout = self.accumulator

        ######### Initialize Objects #########
        jets = events.Jet 
        njets = ak.num(jets)

        ######### Systematics #########
        if not isData: 
            # If this is not an eft sample, get the genWeight
            if eft_coeffs is None: 
                genw = events['genWeight']
            else: 
                genw = np.ones_like(events['event'])

            norm = genw*(xsec/sow)
            weight = norm

        else: 
            weight = np.ones_like(events['event'])

        ######### Add variables to EVENTS #########
        events['njets'] = ak.num(jets)

        ######### Fill dense axes variables ##########
        dense_axis_variables = {}
        dense_axis_variables['njets'] = njets
   
        ########## Fill the histograms ##########
        for dense_axis_name, dense_axis_vals in dense_axis_variables.items():                        

            # Fill the histos
            axes_fill_info_dict = {
                dense_axis_name : dense_axis_vals,
                "process"       : histAxisName,
                "weight"        : weight,
                "eft_coeff"     : eft_coeffs,
            }

            hout[dense_axis_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator


