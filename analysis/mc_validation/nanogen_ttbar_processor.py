#!/usr/bin/env python
import numpy as np
import awkward as ak
import json
import hist
from hist import Hist

from coffea import processor
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

import topcoffea.modules.eft_helper as efth

from ttbarEFT.modules.analysis_tools import get_lumi
from ttbarEFT.modules.gen_tools import is_clean, order_wc_values, calc_event_weights

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


SM_pt = {"ctGIm": 0.0, "ctGRe":0.0, 
         "ctu1": 0.0, "ctd1": 0.0, "ctj1": 0.0, "cQu1": 0.0, "cQd1": 0.0, 
         "ctu8": 0.0, "ctd8": 0.0, "ctj8": 0.0, "cQu8": 0.0, "cQd8": 0.0, 
         "cQj11": 0.0, "cQj18": 0.0,
         "cQj31": 0.0, "cQj38": 0.0,}

EFTpt1 = {"ctGIm": 0.5, "ctGRe":0.5, 
         "ctu1": 1.0, "ctd1": 1.0, "ctj1": 1.0, "cQu1": 1.0, "cQd1": 1.0, 
         "ctu8": 1.0, "ctd8": 1.0, "ctj8": 1.0, "cQu8": 1.0, "cQd8": 1.0, 
         "cQj11": 1.0, "cQj18": 1.0,
         "cQj31": 1.0, "cQj38": 1.0,}

EFTpt2 = {"ctGIm": 2.0, "ctGRe":2.0, 
         "ctu1": 3.0, "ctd1": 3.0, "ctj1": 3.0, "cQu1": 3.0, "cQd1": 3.0, 
         "ctu8": 3.0, "ctd8": 3.0, "ctj8": 3.0, "cQu8": 3.0, "cQd8": 3.0, 
         "cQj11": 3.0, "cQj18": 3.0,
         "cQj31": 3.0, "cQj38": 3.0,}

rwgt_choice = SM_pt
# rwgt_choice = EFTpt1
# rwgt_choice = EFTpt2

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], dtype=np.float32):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 

        print(f"\n\n")
        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)
        print(f"\n\n")

        # Create the histograms with new scikit hist
        self._histo_dict = {
            "njets"         : Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="njets", label='njets'), storage="weight"),
            "ljets"         : Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="ljets", label='ljets'), storage="weight"),
            "nleps"         : Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="nleps", label='nleps'), storage="weight"),
            "ntops"         : Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="ntops", label='ntops'), storage="weight"),
            "mtt"           : Hist(hist.axis.Regular(bins=30, start=0, stop=1500, name='mll', label='mll [GeV]'), storage="weight"),
            "pttt"          : Hist(hist.axis.Regular(bins=35, start=0, stop=700, name='pttt', label='$p_T(tt)$ [GeV]'), storage='weight'),
            "top1pt"        : Hist(hist.axis.Regular(bins=35, start=0, stop=700, name='top1pt', label='top1 $p_T$ [GeV]'), storage='weight'),
            "top2pt"        : Hist(hist.axis.Regular(bins=35, start=0, stop=700, name='top2pt', label='top2 $p_T$ [GeV]'), storage='weight'),
            "top1eta"       : Hist(hist.axis.Regular(bins=30, start=-3, stop=3, name='top1eta', label='top1 eta'), storage='weight'),
            "top2eta"       : Hist(hist.axis.Regular(bins=30, start=-3, stop=3, name='top2eta', label='top2 eta'), storage='weight'),
            "top1phi"       : Hist(hist.axis.Regular(bins=40, start=-4, stop=4, name='top1phi', label='top1 phi'), storage='weight'),
            "top2phi"       : Hist(hist.axis.Regular(bins=40, start=-4, stop=4, name='top2phi', label='top2 phi'), storage='weight'),
            "top1mass"      : Hist(hist.axis.Regular(bins=30, start=140, stop=200, name='top1mass', label='top1 mass [GeV]'), storage='weight'),
            "top2mass"      : Hist(hist.axis.Regular(bins=30, start=140, stop=200, name='top2mass', label='top2 mass [GeV]'), storage='weight'),
            "mll"           : Hist(hist.axis.Regular(bins=40, start=0, stop=800, name='mll', label='mll [GeV]'), storage='weight'),
            "ptll"          : Hist(hist.axis.Regular(bins=25, start=0, stop=500, name='ptll', label='$p_T(ll)$ [GeV]'), storage='weight'),
            "l0pt"          : Hist(hist.axis.Regular(bins=40, start=0, stop=400, name='l0pt', label="l0 $p_T$ [GeV]"), storage="weight"),
            "l1pt"          : Hist(hist.axis.Regular(bins=40, start=0, stop=400, name='10pt', label="l1 $p_T$ [GeV]"), storage="weight"),
            "l0eta"         : Hist(hist.axis.Regular(bins=30, start=-3, stop=3, name='l0eta', label='l0 eta'), storage='weight'),
            "l1eta"         : Hist(hist.axis.Regular(bins=30, start=-3, stop=3, name='l1eta', label='l1 eta'), storage='weight'),
            "l0phi"         : Hist(hist.axis.Regular(bins=40, start=-4, stop=4, name='l0phi', label='l0 phi'), storage='weight'),
            "l1phi"         : Hist(hist.axis.Regular(bins=40, start=-4, stop=4, name='l1phi', label='l1 phi'), storage='weight'),
            "dr_leps"       : Hist(hist.axis.Regular(bins=24, start=0, stop=6, name='dr_leps', label="Delta R (leading, subleading lepton)"), storage="weight"),
            "j0pt"          : Hist(hist.axis.Regular(bins=50, start=0, stop=500, name='j0pt', label='j0$p_T$ [GeV]'), storage='weight'),
            "j0eta"         : Hist(hist.axis.Regular(bins=30, start=-3, stop=3, name='j0eta', label='j0 eta'), storage='weight'),
            "j0phi"         : Hist(hist.axis.Regular(bins=40, start=-4, stop=4, name='j0phi', label='j0 phi'), storage='weight'),
        }

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


        # Initialize Objects
        genpart = events.GenPart
        is_final_mask = genpart.hasFlags(["fromHardProcess","isLastCopy"])
        nu_ele = genpart[is_final_mask & (abs(genpart.pdgId) == 12)]
        nu_mu = genpart[is_final_mask & (abs(genpart.pdgId) == 14)]
        ele  = genpart[is_final_mask & (abs(genpart.pdgId) == 11)]
        mu   = genpart[is_final_mask & (abs(genpart.pdgId) == 13)]


        ######## Lepton Selections ########
        nu = ak.concatenate([nu_ele,nu_mu],axis=1)
        e_selec = ((ele.pt>20) & (abs(ele.eta)<2.5))
        m_selec = ((mu.pt>20) & (abs(mu.eta)<2.5))

        # leps = ak.concatenate([ele, mu], axis=1)
        leps = ak.concatenate([ele[e_selec],mu[m_selec]],axis=1)
        leps = leps[ak.argsort(leps.pt, axis=-1, ascending=False)]


        ######## Jet Selections ########
        jets = events.GenJet
        jets = jets[(jets.pt>30) & (abs(jets.eta)<2.5)]
        jets_clean = jets[is_clean(jets, leps, drmin=0.4) & is_clean(jets, nu, drmin=0.4)]
        j0 = jets_clean[ak.argmax(jets_clean.pt, axis=-1, keepdims=True)]


        ######## Top Selections ########
        gen_top = ak.pad_none(genpart[is_final_mask & (abs(genpart.pdgId) == 6)],2)
        gen_top = gen_top[ak.argsort(gen_top.pt, axis=1, ascending=False)]
        mtt = (gen_top[:,0] + gen_top[:,1]).mass


        ######## Event selections ########
        nleps = ak.num(leps)
        ntops = ak.num(gen_top)
        njets = ak.num(jets_clean)
        ljets = ak.num(jets_clean[abs(jets_clean.partonFlavour) != 5])
        
        at_least_two_leps = ak.fill_none(nleps>=2,False)
        at_least_two_jets = ak.fill_none(njets>=2,False)

        selections = PackedSelection()
        selections.add('2l', at_least_two_leps)
        selections.add('2j', at_least_two_jets)

        # event_selection_mask = selections.all('2l')
        event_selection_mask = selections.all('2l', '2j')


        ######## Variables to Plot ########
        leps = ak.pad_none(leps, 2)
        l0 = leps[:,0]
        l1 = leps[:,1]

        ptll = (l0+l1).pt
        dr_leps = l0.delta_r(l1)

        # mtt = (gen_top[:,0]+gen_top[:,1]).mass
        # pttt = (gen_top[:,0]+gen_top[:,1]).pt


        ######## Normalization ########
        lumi = 1000.0*get_lumi(year) #lumi is in fb^-1 but xsec is in pb
        norm = (xsec/sow)*lumi 

        counts = np.ones_like(events['event'])[event_selection_mask]

        # eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

        if eft_coeffs is None: 
            event_weights = events["genWeight"]
        else:  
            wc_lst = order_wc_values(self._wc_names_lst, rwgt_choice)
            event_weights = calc_event_weights(eft_coeffs, wc_lst)


        ######## Fill Histograms ########
        hout = self._histo_dict

        variables_to_fill = {}
        variables_to_fill['njets'] = njets
        variables_to_fill['ljets'] = ljets
        # variables_to_fill['nleps'] = nleps
        # variables_to_fill['ntops'] = ntops
        
        variables_to_fill['mtt'] = (gen_top[:,0]+gen_top[:,1]).mass
        variables_to_fill['pttt'] = (gen_top[:,0]+gen_top[:,1]).pt
        variables_to_fill['top1pt'] = gen_top[:,0].pt 
        variables_to_fill['top2pt'] = gen_top[:,1].pt
        variables_to_fill['top1eta'] = gen_top[:,0].eta
        variables_to_fill['top2eta'] = gen_top[:,1].eta
        variables_to_fill['top1phi'] = gen_top[:,0].phi
        variables_to_fill['top2phi'] = gen_top[:,1].phi
        variables_to_fill['top1mass'] = gen_top[:,0].mass
        variables_to_fill['top2mass'] = gen_top[:,1].mass

        variables_to_fill['mll'] = (l0+l1).mass
        variables_to_fill['ptll'] = (l0+l1).pt
        variables_to_fill['l0pt'] = l0.pt
        variables_to_fill['l1pt'] = l1.pt
        variables_to_fill['l0eta'] = l0.eta
        variables_to_fill['l1eta'] = l1.eta
        variables_to_fill['l0phi'] = l0.phi
        variables_to_fill['l1phi'] = l1.phi

        variables_to_fill['dr_leps'] = dr_leps
        variables_to_fill['j0pt'] = ak.flatten(j0.pt) 
        variables_to_fill['j0eta'] = ak.flatten(j0.eta)
        variables_to_fill['j0phi'] = ak.flatten(j0.phi)

        for var_name, var_values in variables_to_fill.items():
            # if var_name not in self._hist_lst:
            #     print(f"Skipping \"{var_name}\", it is not in the list of hists to include")
            #     continue

            hout[var_name].fill(var_values[event_selection_mask], weight=(event_weights[event_selection_mask])*norm)

        return hout

    def postprocess(self, accumulator):
        return accumulator        

