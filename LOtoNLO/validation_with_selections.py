import awkward as ak
import numpy as np
import pandas as pd


from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.analysis_tools import PackedSelection

import hist
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth

import ttbarEFT.modules.object_selection as tt_os
from ttbarEFT.modules.analysis_tools import make_mt2, get_lumi

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


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

        axes = {
            "reweights": {
                "regular": (200, 0, 100),
                "label": "reweight values",},
            "sow": {
                "regular": (1, 0, 2),
                "label": "sum of weights",},
            "ptll": {
                "regular": (25, 0, 500),
                "label": "$p_T(ll)$ [GeV]"},
            "pttt": {
                "regular": (35, 0, 700),
                "label": "$p_T(tt)$ [GeV]"},
            "mtt": {
                "regular": (75, 0, 1500),
                "label": "mtt"},
            "lep1pt": {
                "regular": (40, 0, 400),
                "label": "lep1 pt",},
            "lep2pt": {
                "regular": (40, 0, 400),
                "label": "lep2 pt",},
            "top1pt": {
                "regular": (35, 0, 700),
                "label": "top1 pt",},
            "top2pt": {
                "regular": (35, 0, 700),
                "label": "top2 pt",},
            "lep1eta": {
                "regular": (50, -5, 5),
                "label": "lep1 eta",},
            "lep2eta": {
                "regular": (50, -5, 5),
                "label": "lep2 eta",},
            "top1eta": {
                "regular": (50, -5, 5),
                "label": "top1 eta",},
            "top2eta": {
                "regular": (50, -5, 5),
                "label": "top2 eta",},
            "lep1phi": {
                "regular": (40, -4, 4),
                "label": "lep1 phi",},
            "lep2phi": {
                "regular": (40, -4, 4),
                "label": "lep2 phi",},
            "top1phi": {
                "regular": (40, -4, 4),
                "label": "top1 phi",},
            "top2phi": {
                "regular": (40, -4, 4),
                "label": "top2 phi",},
            "top1mass": {
                "regular": (34, 80, 250),
                "label": "top1 mass", },
            "top2mass": {
                "regular": (34, 80, 250),
                "label": "top2 mass", },
            "j0pt": {
                "regular": (100, 0, 500),
                "label": "j0pt",},
            "j0eta": {
                "regular": (50, -5, 5),
                "label": "j0eta",},
            "j0phi": {
                "regular": (40, -4, 4),
                "label": "j0phi",},
            "njets": {
                "regular": (10, 0, 10),
                "label": "njets",},
        }

        histograms = {}

        for name, info in axes.items():
            if 'variable' in info:
                dense_axis = hist.axis.Variable(info['variable'], name=name, label=info['label'])
            else:
                dense_axis = hist.axis.Regular(*info['regular'], name=name, label=info['label'])

            histograms[name]=HistEFT(
                proc_axis,
                dense_axis,
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


        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None
        genpart = events.GenPart
        is_final_mask = genpart.hasFlags(["fromHardProcess","isLastCopy"])

        gen_top = ak.pad_none(genpart[is_final_mask & (abs(genpart.pdgId) == 6)],2)
        gen_top = gen_top[ak.argsort(gen_top.pt, axis=1, ascending=False)]

        ele  = genpart[is_final_mask & (abs(genpart.pdgId) == 11)]
        mu   = genpart[is_final_mask & (abs(genpart.pdgId) == 13)]
        tau  = genpart[is_final_mask & (abs(genpart.pdgId) == 15)]
        nu_ele = genpart[is_final_mask & (abs(genpart.pdgId) == 12)]
        nu_mu = genpart[is_final_mask & (abs(genpart.pdgId) == 14)]
        nu_tau = genpart[is_final_mask & (abs(genpart.pdgId) == 16)]
        nu = ak.concatenate([nu_ele,nu_mu, nu_tau],axis=1)

        e_selec = ((ele.pt>20) & (abs(ele.eta)<2.5))
        m_selec = ((mu.pt>20) & (abs(mu.eta)<2.5))
        t_selec = ((tau.pt>20) & (abs(tau.eta)< 2.5))
        leps = ak.concatenate([ele[e_selec], mu[m_selec], tau[t_selec]],axis=1)

        leps = leps[ak.argsort(leps.pt, axis=-1, ascending=False)]
        nleps = ak.num(leps)

        jets = events.GenJet
        jets = jets[(jets.pt>30) & (abs(jets.eta)<2.5)]
        jets_clean = jets[tt_os.isClean(jets, leps, drmin=0.4) & tt_os.isClean(jets, nu, drmin=0.4)]
        jets_clean = jets_clean[ak.argsort(jets_clean.pt, axis=-1, ascending=False)]

        jets = jets_clean[ak.argsort(jets_clean.pt, axis=-1, ascending=False)]
        j0 = jets_clean[ak.argmax(jets_clean.pt, axis=-1, keepdims=True)]
        njets = ak.num(jets_clean)
        padded_jets = ak.pad_none(jets_clean, 3, axis=1)

        ######## Event selections ########

        selections = PackedSelection()

        at_least_two_leps = ak.fill_none(nleps>=2, False)
        at_least_two_jets = ak.fill_none(njets>=2, False)

        selections.add('2l', at_least_two_leps)
        selections.add('2j', at_least_two_jets)

        event_selection_mask = selections.all('2l', '2j')

        ######## Variables for Plotting ########

        leps = ak.pad_none(leps, 2)
        l0 = leps[:,0]
        l1 = leps[:,1]

        ptll = (l0+l1).pt

        ######## Normalizations ########

        lumi = 1000.0*get_lumi(year)

        norm = (xsec/sow)*lumi

        if eft_coeffs is None:
            genw = events["genWeight"]
        else:
            genw = np.ones_like(events['event'])

        event_weights = norm*genw

        # NLO_weight = tt_cor.GetNLO_Weight(events, dataset)
        # event_weights = norm*genw*NLO_weight

        ######## Fill Histograms ########

        hout = self.accumulator

        variables_to_fill = {
            "sow"       : np.ones_like(events['event']),
            "pttt"      : (gen_top[:,0] + gen_top[:,1]).pt,
            "mtt"       : (gen_top[:,0] + gen_top[:,1]).mass,
            "top1pt"    : gen_top.pt[:,0],
            "top1eta"   : gen_top.eta[:,0],
            "top1phi"   : gen_top.phi[:,0],
            "top1mass"  : gen_top.mass[:,0],
            "top2pt"    : gen_top.pt[:,1],
            "top2eta"   : gen_top.eta[:,1],
            "top2phi"   : gen_top.phi[:,1],
            "top2mass"  : gen_top.mass[:,1],
            "ptll"      : ptll,
            "lep1pt"    : l0.pt,
            "lep1eta"   : l0.eta,
            "lep1phi"   : l0.phi,
            "lep2pt"    : l1.pt,
            "lep2eta"   : l1.eta,
            "lep2phi"   : l1.phi,
            "j0pt"      : ak.fill_none(padded_jets.pt[:,0], 0.0),
            "j0eta"     : ak.fill_none(padded_jets.eta[:,0], 0.0),
            "j0phi"     : ak.fill_none(padded_jets.phi[:,0], 0.0),
            "njets"     : njets,
        }     

        eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

        for var_name, var_values in variables_to_fill.items():

            fill_info = {
                var_name    : var_values[event_selection_mask],
                "process"   : histAxisName,
                "weight"    : event_weights[event_selection_mask],
                "eft_coeff" : eft_coeffs_cut,
            }

            hout[var_name].fill(**fill_info)

        return hout

    def postprocess(self, accumulator):
        return accumulator

