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

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth
import topcoffea.modules.corrections as tc_cor
import topcoffea.modules.event_selection as tc_es

from ttbarEFT.modules.analysis_tools import make_mt2, get_lumi
from ttbarEFT.modules.axes import info as axes_info
import ttbarEFT.modules.event_selection as tt_es
import ttbarEFT.modules.object_selection as tt_os

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systs=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, muonSyst='nominal', dtype=np.float32, rebin=False, offZ_split=False, tau_h_analysis=False, fwd_analysis=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)
        appl_axis = hist.axis.StrCategory([], name="appl", label=r"AR/SR", growth=True)

        # Set the booleans
        self._do_errors = do_errors # Whether to calculate and store the w**2 coefficients
        self._do_systematics = do_systs # Whether to process systematic samples
        self._split_by_lepton_flavor = split_by_lepton_flavor # Whether to keep track of lepton flavors individually
        self._skip_signal_regions = skip_signal_regions # Whether to skip the SR categories
        self._skip_control_regions = skip_control_regions # Whether to skip the CR categories

        # fill histograms using info from axes.py 
        histograms = {}
        for name, info in axes_info.items():
            if "variable" in info: 
                dense_axis = hist.axis.Variable(info["variable"], name=name, label=info["label"])
                sumw2_axis = hist.axis.Variable(info["variable"], name=name+"_sumw2", label=info["label"]+" sum of w^2")
            else:
                dense_axis = hist.axis.Regular(*info["regular"], name=name, label=info["label"])
                sum2w_axis = hist.axis.Regular(*info["regular"], name=name+"_sumw2", label=info["label"] + " sum of w^2")

            histograms[name] = HistEFT(
                proc_axis, 
                syst_axis,
                dense_axis,
                wc_names = wc_names_lst, 
                label=r"Events",
            )

            histograms[name+"_sumw2"] = HistEFT(
                proc_axis, 
                syst_axis,
                sum2w_axis,
                wc_names = wc_names_lst, 
                label=r"Events",
            )

        self._accumulator = histograms

        # set the list of hists to fill
        if hist_lst is None:
            self._hist_lst = list(self._accumulator.keys()) #fill all hists if not specified
        else:
            for hist_to_include in hist_lst:
                if hist_to_include not in self._accumulator.keys():
                    raise Exception(f"Error: Cannot specify hist \"{hist_to_include}\", it is not defined in the processor.")
            self._hist_lst = hist_lst 

        # print out basic info before running over files
        print("\n\n")
        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)
        print("hist_lst: ", self._hist_lst)
        print("\n\n")

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns


    def process(self, events):

        # Dataset parameters
        dataset         = events.metadata['dataset']
        # isEFT             = self._samples[dataset]["WCnames"] != []
        isData          = self._samples[dataset]["isData"]
        histAxisName  = self._samples[dataset]["histAxisName"]
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']

        isEFT = hasattr(events, "EFTfitCoefficients")    
        assert not (isEFT and isData), f"isEFT and isData cannot both be True. Check input samples."

        
        ######### Initialize Objects #########

        met  = events.MET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet 

        # An array of lenght events that is just 1 for each event
        events.nom = ak.ones_like(events.MET.pt)

        leptonSelection = tt_os.Run2LeptonSelection()

        ######### Leptons  ##########
        ele['isPres']=leptonSelection.is_pres_ele(ele)
        mu['isPres']=leptonSelection.is_pres_muon(mu)

        ele_good = ele[ele.isPres]
        mu_good = mu[mu.isPres]

        leps = ak.concatenate([ele_good, mu_good], axis=1)
        nleps = ak.num(leps)

        # leps = ak.concatenate([ele,mu],axis=1)
        # nleps = ak.num(leps)


        ######### EFT coefficients ##########

        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events["EFTfitCoefficients"]) if hasattr(events, "EFTfitCoefficients") else None
        if eft_coeffs is not None:
            # Check to see if the ordering of WCs for this sample matches what want
            if self._samples[dataset]["WCnames"] != self._wc_names_lst:
                eft_coeffs = efth.remap_coeffs(self._samples[dataset]["WCnames"], self._wc_names_lst, eft_coeffs)
        eft_w2_coeffs = efth.calc_w2_coeffs(eft_coeffs,self._dtype) if (self._do_errors and eft_coeffs is not None) else None
        # Initialize the out object
        hout = self.accumulator


        ######### Selection Masks #########

        # create trigger mask
        pass_trg = tc_es.trg_pass_no_overlap(events, isData, dataset, str(year), tt_es.triggers_dict, tt_es.exclude_triggers_dict)


        ######### Store boolean masks with PackedSelection ##########

        selections = PackedSelection(dtype='uint64')
        
        at_least_two_leps = ak.fill_none(nleps>=2, False)
        selections.add("2l", at_least_two_leps) 
        selections.add("trg", pass_trg)

        event_selection_mask = selections.all('2l', 'trg')

        good_events = events[event_selection_mask]

        fname = 'output_lep_selec'
        with open(f'{fname}.txt', 'w') as output:
            for i in range(len(good_events)):
                output.write(f"{good_events[i].run}:{get_lumi(year)}:{good_events[i].event}\n")

        print(f"\n\n run:lumi:event info saved to {fname}.txt\n\n ")

        return hout

    def postprocess(self, accumulator):
        return accumulator


