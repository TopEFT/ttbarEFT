#!/usr/bin/env python
import copy
import coffea
import numpy as np
import awkward as ak
import json
import hist
# import yaml

from coffea import processor
from coffea.analysis_tools import PackedSelection
# from coffea.nanoevents.methods import vector
# from coffea.lumi_tools import LumiMask

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

# from topcoffea.modules.paths import topcoffea_path
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth
# import topcoffea.modules.corrections as tc_cor
# import topcoffea.modules.event_selection as tc_es

from ttbarEFT.modules.paths import ttbarEFT_path
from ttbarEFT.modules.analysis_tools import get_lumi
# import ttbarEFT.modules.object_selection as tt_os
# import ttbarEFT.modules.event_selection as tt_es
# from ttbarEFT.modules.corrections import AttachElectronSF, AttachMuonSF, ApplyMuonPtCorr, ApplyJetVetoMaps, AttachElecTrigEff, AttachMuonTrigEff, AttachTrigSF
# import ttbarEFT.modules.corrections as tt_cor 

# from topcoffea.modules.get_param_from_jsons import GetParam

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, lep_cat, wc_names_lst=[], hist_lst=None, dtype=np.float32):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 
        self._lep_cat = lep_cat

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)

        histograms = {}
        histograms['mtt'] = HistEFT(
                                proc_axis,
                                hist.axis.Regular(bins=125, start=250, stop=1500, name='mtt', label='mtt'),
                                wc_names = wc_names_lst,
                                label=r'Events')


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

        lhepart = events.LHEPart
        lhe_mtt = (lhepart[:,0]+lhepart[:,1]+lhepart[:,2]+lhepart[:,3]+lhepart[:,4]+lhepart[:,5]+lhepart[:,6]+lhepart[:,7]).mass


        if eft_coeffs is None:
            genw = events["genWeight"]
        else:
            genw = np.ones_like(events['event'])

        lumi = 1000.0*get_lumi(year)
        norm = (xsec/sow)*lumi

        variables_to_fill = {
            "mtt"   : lhe_mtt,
        }

        for var_name, var_values in variables_to_fill.items():
            fill_info = {
                var_name    : var_values,
                "process"   : histAxisName,
                "weight"    : norm*genw,
                "eft_coeff" : eft_coeffs,
            }

            hout[var_name].fill(**fill_info)

        return hout

    def postprocess(self, accumulator):
        return accumulator

