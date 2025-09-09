#!/usr/bin/env python
import copy
import coffea
import numpy as np
import awkward as ak
import json
import hist

from mt2 import mt2

from coffea import processor
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods import vector
from coffea.lumi_tools import LumiMask

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from topcoffea.modules.paths import topcoffea_path
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth
import topcoffea.modules.corrections as tc_cor
import topcoffea.modules.event_selection as tc_es

from ttbarEFT.modules.paths import ttbarEFT_path
from ttbarEFT.modules.analysis_tools import make_mt2, get_lumi
import ttbarEFT.modules.object_selection as tt_os
import ttbarEFT.modules.event_selection as tt_es
from ttbarEFT.modules.corrections import AttachElectronSF, AttachMuonSF, AttachMuonTrigSF, ApplyMuonPtCorr, ApplyJetVetoMaps

from topcoffea.modules.get_param_from_jsons import GetParam

get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_tt_param = GetParam(ttbarEFT_path("params/params.json"))

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], dtype=np.float32):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)

        # fill histograms using info from axes.json
        with open(ttbarEFT_path("params/axes.json"), 'r') as axes_file:
            axes_info = json.load(axes_file)

        histograms = {}
        for name, info in axes_info['CR_axes'].items():
        # for name, info in axes_info.items():
            if 'variable' in info: 
                dense_axis = hist.axis.Variable(info['variable'], name=name, label=info['label'])
                sumw2_axis = hist.axis.Variable(info['variable'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')
            else:
                dense_axis = hist.axis.Regular(*info['regular'], name=name, label=info['label'])
                sumw2_axis = hist.axis.Regular(*info['regular'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')

            histograms[name] = HistEFT(
                proc_axis, 
                # chan_axis,
                # syst_axis,
                dense_axis,
                wc_names = wc_names_lst, 
                label=r'Events',
            )

            histograms[name+'_sumw2'] = HistEFT(
                proc_axis, 
                # chan_axis,
                # syst_axis,
                sumw2_axis,
                wc_names = wc_names_lst, 
                label=r'Events',
            )

        self._accumulator = histograms

        # print out basic info before running over filesrun2leptonselection
        print("\n\n")
        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)
        print("\n\n")

    @property
    def accumulator(self):
        return self._accumulator

    # @property
    # def columns(self):
    #     return self._columns

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

        ######### Lumi Mask for Data #########        
        golden_json_path = {
            "2016": topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
            "2016APV": topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
            "2017": topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
            "2018": topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
        }
        lumi_mask = LumiMask(golden_json_path[year])(events.run,events.luminosityBlock)


        ######### Initialize Objects #########
        met  = events.MET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet 

        leptonSelection = tt_os.Run2LeptonSelection()

        # An array of length events that is just 1 for each event
        events['nom'] = ak.ones_like(events.MET.pt)

        ######### Electron Selection ##########
        ele['isSelE']=leptonSelection.is_sel_ele(ele)
        ele_good = ele[ele.isSelE]

        AttachElectronSF(ele_good, year)    # TODO: this should be added inside an only MC block


        ######### Muon Selection ##########
        mu['pt'] = ApplyMuonPtCorr(mu, year, isData)
        mu['isSelM']=leptonSelection.is_sel_muon(mu)
        
        mu_good = mu[mu.isSelM]
        AttachMuonSF(mu_good, year)         # TODO: this should be added inside an only MC block


        ######### Lepton Selection ##########
        leps = ak.concatenate([ele_good, mu_good], axis=1)
        leps_sorted = leps[ak.argsort(leps.pt, axis=-1,ascending=False)] 

        events['leps_pt_sorted'] = leps_sorted


        ######### Add variables into event object so that they persist #########
        events['njets'] = ak.num(jets)

        ######### Fill dense axes variables ##########
        dense_axis_variables = {}
        dense_axis_variables['njets'] = events.njets
        weights = np.ones_like(events['event'])


        ########## Fill the histograms ##########

        for dense_axis_name, dense_axis_vals in dense_axis_variables.items():

            # Fill the histos
            axes_fill_info_dict = {
                dense_axis_name : dense_axis_vals,
                "process"       : histAxisName,
                "weight"        : weights,
                "eft_coeff"     : eft_coeffs,
            }

            hout[dense_axis_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator