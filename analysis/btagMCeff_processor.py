#!/usr/bin/env python
import copy
import coffea
import numpy as np
import awkward as ak
import json
import hist
import yaml

from coffea import processor
from coffea.analysis_tools import PackedSelection

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from topcoffea.modules.paths import topcoffea_path
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth
# import topcoffea.modules.corrections as tc_cor
# import topcoffea.modules.event_selection as tc_es

from ttbarEFT.modules.paths import ttbarEFT_path
import ttbarEFT.modules.object_selection as tt_os
import ttbarEFT.modules.event_selection as tt_es
import ttbarEFT.modules.corrections as tt_cor 

from topcoffea.modules.get_param_from_jsons import GetParam
get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_tt_param = GetParam(ttbarEFT_path("params/params.json"))

# silence warnings due to using NanoGEN instead of full NanoAOD
NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, lep_cat, wc_names_lst=[], hist_lst=None, dtype=np.float32):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 
        self._lep_cat = lep_cat

        proc_axis = hist.axis.StrCategory([], name='process', growth=True)
        jet_pt_axis = hist.axis.Variable([20, 30, 50, 70, 100, 140, 200, 300, 600, 1000], name='jpt', label=r'Jet p_{T} [GeV]')
        jet_eta_axis = hist.axis.Variable([0, 0.6, 1.2, 2.4], name='jeta', label=r'Jet \eta [GeV]')
        jet_flav_axis = hist.axis.StrCategory([], name='flav', growth=True)
        # jet_flavor_axis = hist.axis.Regular(5, 0, 5, name="flavour", label="jet flavour (int)")
        # jet_flavor_axis = hist.axis.IntCategory([0,1,2,3,4,5], name="flavour", label="jet flavour (int)")
        jet_flavor_axis = hist.axis.IntCategory([0,4,5], name="flavour", label="jet flavour (int)")
        wp_axis = hist.axis.StrCategory([], name="WP", growth=True)

        histograms = {}
        histograms['jetpt'] = hist.Hist(wp_axis, proc_axis, jet_flav_axis, jet_pt_axis)
        histograms['jeteta'] = hist.Hist(wp_axis, proc_axis, jet_flav_axis, jet_eta_axis)
        histograms['jetpteta'] = hist.Hist(wp_axis, proc_axis, jet_flav_axis, jet_pt_axis, jet_eta_axis)
        histograms['jetptetaflav'] = hist.Hist(wp_axis, jet_flavor_axis, jet_pt_axis, jet_eta_axis)

        self._accumulator = histograms

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    def process(self, events):
        dataset         = events.metadata['dataset']
        isData          = self._samples[dataset]['isData']
        histAxisName    = self._samples[dataset]['histAxisName']
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']

        datasets = ["Muon", "SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG"]
        for d in datasets:
            if dataset.startswith(d):
                dataset = dataset.split('_')[0]

        hout = self.accumulator


        ######### Initialize Objects #########
        met  = events.MET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet 

        ######### Object Selection ######### 

        leptonSelection = tt_os.Run2LeptonSelection()

        ele['isGoodElec']=leptonSelection.is_sel_ele(ele)
        ele_good = ele[ele.isGoodElec]
        mu['pt'] = tt_cor.ApplyMuonPtCorr(mu, year, isData)
        mu['isGoodMuon']=leptonSelection.is_sel_muon(mu)
        mu_good = mu[mu.isGoodMuon]

        leps = ak.concatenate([ele_good, mu_good], axis=1)
        leps_sorted = leps[ak.argsort(leps.pt, axis=-1,ascending=False)] 
        events['leps_pt_sorted'] = leps_sorted

        jets['isPres'] = tt_os.is_pres_jet(jets) 
        jets['isClean'] = tt_os.isClean(jets, ele_good, drmin=0.4)& tt_os.isClean(jets, mu_good, drmin=0.4)
        # jets['isClean'] = tt_os.isClean(jets, leps_sorted, drmin=0.4)
        goodJets = jets[(jets.isPres)&jets.isClean] 
        goodJets = goodJets[(goodJets.partonFlavour != 0)]

        medium_tag = "btag_wp_medium_" + year.replace("201", "UL1")
        btagwpm = get_tc_param(medium_tag)
        isBtagJetsMedium = (goodJets.btagDeepFlavB > btagwpm)

        tt_es.add2losMask(events, year, isData)
        selections = PackedSelection(dtype='uint64')
        selections.add('2los', events.is2los)

        btagSelection = {
            'all' : goodJets.btagDeepFlavB > -999.0, 
            'medium': goodJets.btagDeepFlavB > btagwpm,
        }

        flavSelection = {
            'b': (np.abs(goodJets.hadronFlavour) == 5), 
            'c': (np.abs(goodJets.hadronFlavour) == 4), 
            'l': (np.abs(goodJets.hadronFlavour) <= 3), 
        }

        # for jetflav in flavSelection.keys():
        #     for wp in btagSelection.keys():
        #         mask = (flavSelection[jetflav])&(btagSelection[wp])&(events.is2los)
        #         selectedJets = goodJets[mask]

        #         pts = ak.flatten(selectedJets.pt)
        #         etas = ak.flatten(selectedJets.eta)
        #         absetas = ak.flatten(np.abs(selectedJets.eta))

        #         flavarray = np.zeros_like(pts) if jetflav == 'l' else (np.ones_like(pts)*(4 if jetflav=='c' else 5))
        #         weights = np.ones_like(pts)

        #         hout['jetpt'].fill(WP=wp, process=histAxisName, flav=jetflav, jpt=pts, weight=weights)
        #         hout['jeteta'].fill(WP=wp, process=histAxisName, flav=jetflav, jeta=etas, weight=weights)
        #         hout['jetpteta'].fill(WP=wp, process=histAxisName, flav=jetflav, jpt=pts, jeta=absetas, weight=weights)
        #         hout['jetptetaflav'].fill(WP=wp, process=histAxisName, flavour=flavarray, jpt=pts, jeta=absetas, weight=weights)

        for wp in btagSelection.keys():
            mask = (btagSelection[wp])&(events.is2los)
            selectedJets = goodJets[mask]

            pts = ak.flatten(selectedJets.pt)
            etas = ak.flatten(selectedJets.eta)
            absetas = ak.flatten(np.abs(selectedJets.eta))
            flavarray = ak.flatten(selectedJets.hadronFlavour)

            # flavarray = np.zeros_like(pts) if jetflav == 'l' else (np.ones_like(pts)*(4 if jetflav=='c' else 5))
            # weights = np.ones_like(pts)

            # hout['jetpt'].fill(WP=wp, process=histAxisName, flav=jetflav, jpt=pts, weight=weights)
            # hout['jeteta'].fill(WP=wp, process=histAxisName, flav=jetflav, jeta=etas, weight=weights)
            # hout['jetpteta'].fill(WP=wp, process=histAxisName, flav=jetflav, jpt=pts, jeta=absetas, weight=weights)
            hout['jetptetaflav'].fill(WP=wp, jpt=pts, jeta=absetas, flavour=flavarray)


        return hout


    def postprocess(self, accumulator):
        return accumulator

