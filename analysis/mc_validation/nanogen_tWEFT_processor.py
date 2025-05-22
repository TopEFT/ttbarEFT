#!/usr/bin/env python
import numpy as np
import awkward as ak

np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods import vector

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = False

import hist
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth

from mt2 import mt2
# from mt2 import mt2_arxiv

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

def make_mt2(l0, l1, met):
    nevents = len(np.zeros_like(met))
    misspart = ak.zip(
        {
            "pt": met.pt,
            "eta": 0,
            "phi": met.phi,
            "mass": np.full(nevents, 0),
        },
        with_name="PtEtaPhiMLorentzVector",
        behavior=vector.behavior,
    )

    # mt2_val = mt2_arxiv(...)
    mt2_val = mt2(
        l0.mass, l0.px, l0.py,                          # visible particle #1
        l1.mass, l1.px, l1.py,                          # visible particle #2 
        misspart.px, misspart.py,                       # missing transverse momentum
        np.zeros_like(met.pt), np.zeros_like(met.pt)    # invisible masses
    )

    return mt2_val

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
            "njets"     :HistEFT(
                            proc_axis, 
                            hist.axis.Regular(bins=10, start=0, stop=10, name="njets", label="njets"),
                            wc_names=wc_names_lst,
                            label="Events"),
            "nleps"     :HistEFT(
                            proc_axis, 
                            hist.axis.Regular(bins=10, start=0, stop=10, name="nleps", label="nleps"),
                            wc_names=wc_names_lst,
                            label="Events"),
            "ntops"     :HistEFT(
                            proc_axis, 
                            hist.axis.Regular(bins=10, start=0, stop=10, name="ntops", label="ntops"),
                            wc_names=wc_names_lst,
                            label="Events"),
            "mll"       :HistEFT(
                            proc_axis, 
                            hist.axis.Regular(bins=40, start=0, stop=800, name='mll', label='invariant mass of the leptons [GeV]'),
                            wc_names=wc_names_lst,
                            label="Events"),
            "dr_leps"   :HistEFT(
                            proc_axis, 
                            hist.axis.Regular(bins=24, start=0, stop=8, name='dr_leps', label='$\Delta R$ (leading lepton, subleading lepton)'),
                            wc_names=wc_names_lst,
                            label="Events"),
            "l0pt"      :HistEFT(
                            proc_axis, 
                            hist.axis.Regular(bins=40, start=0, stop=400, name='l0pt', label='leading lepton $p_T$ [GeV]'),
                            wc_names=wc_names_lst,
                            label="Events"),
            "top_pt"   :HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=70, start=0, stop=700, name="top_pt", label="top $p_T$ [GeV]"),
                            wc_names=wc_names_lst,
                            label="Events"),
            "sow"       :HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow", label="sum of weights for all events"), 
                            wc_names=wc_names_lst, 
                            label="Events"),
            "mt2"       :HistEFT(
                            proc_axis, 
                            hist.axis.Variable([0,10,50,80,100,110,120,140,160,180,220,500], name="mt2", label="$m_{T2}$ [GeV]"),
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
        genpart = events.GenPart
        is_final_mask = genpart.hasFlags(["fromHardProcess","isLastCopy"])
        ele  = genpart[is_final_mask & (abs(genpart.pdgId) == 11)]
        mu   = genpart[is_final_mask & (abs(genpart.pdgId) == 13)]
        nu_ele = genpart[is_final_mask & (abs(genpart.pdgId) == 12)]
        nu_mu = genpart[is_final_mask & (abs(genpart.pdgId) == 14)]
        nu = ak.concatenate([nu_ele,nu_mu],axis=1)
        jets = events.GenJet

        ######## Lep selection  ########

        e_selec = ((ele.pt>20) & (abs(ele.eta)<2.5))
        m_selec = ((mu.pt>20) & (abs(mu.eta)<2.5))
        leps = ak.concatenate([ele[e_selec],mu[m_selec]],axis=1)
        leps = leps[ak.argsort(leps.pt, axis=-1, ascending=False)]
        l0 = leps[ak.argmax(leps.pt, axis=-1, keepdims=True)]

        ######## Jet selection  ########

        jets = jets[(jets.pt>30) & (abs(jets.eta)<2.5)]
        jets_clean = jets[is_clean(jets, leps, drmin=0.4) & is_clean(jets, nu, drmin=0.4)]
        ht = ak.sum(jets_clean.pt, axis=-1)
        j0 = jets_clean[ak.argmax(jets_clean.pt, axis=-1, keepdims=True)]

        ######## Top selection ########

        gen_top = genpart[is_final_mask & (abs(genpart.pdgId) == 6)]

        ######## Event selections ########

        nleps = ak.num(leps)
        njets = ak.num(jets_clean)
        ntops = ak.num(gen_top)

        # Standard 2l2j selections
        at_least_two_leps = ak.fill_none(nleps>=2,False)
        at_least_one_jet = ak.fill_none(njets>=1,False)

        selections = PackedSelection()
        selections.add('2l', at_least_two_leps)
        selections.add('1j', at_least_one_jet)
        event_selection_mask = selections.all('2l', '1j')

        leps_cut = leps[event_selection_mask]

        ######## mt2 variable ########
        met  = events.GenMET[event_selection_mask]
        # met  = events.PuppiMET #for reco later
        mt2_var = make_mt2(leps_cut[:,0], leps_cut[:,1], met)

        ######## Variables with Selections ########

        top = gen_top[event_selection_mask]
        top_pt_cut = top.pt
        
        dr_l0 = leps_cut[:,0]
        dr_l1 = leps_cut[:,1]
        l0pt_cut = l0.pt[event_selection_mask]
        dr_cut = dr_l0.delta_r(dr_l1)

        njets_cut = njets[event_selection_mask]
        nleps_cut = nleps[event_selection_mask]
        ht_cut = ht[event_selection_mask]
        ntops_cut = ntops[event_selection_mask]
        jets_pt_cut = jets_clean.sum().pt[event_selection_mask]
        j0pt_cut = j0.pt[event_selection_mask]
        mll = (leps_cut[:,0] + leps_cut[:,1]).mass

        ######## Normalization ########

        # Normalize by (xsec/sow)
        lumi = 1000.0*get_lumi(year)
        norm = (xsec/sow)*lumi

        if eft_coeffs is None:
            genw = events["genWeight"]
        else:
            genw = np.ones_like(events['event'])

        event_weights = norm*genw

        ######## Fill histos ########

        hout = self._histo_dict

        variables_to_fill = {
            "njets"     : njets_cut,
            "nleps"     : nleps_cut,
            "ntops"     : ntops_cut,
            "mll"       : mll,
            "dr_leps"   : dr_cut,
            "l0pt"      : ak.flatten(l0pt_cut),
            "top_pt"    : ak.flatten(top_pt_cut),
            "mt2"       : mt2_var,
        }

        eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

        for var_name, var_values in variables_to_fill.items():
            if var_name not in self._hist_lst:
                print(f"Skipping \"{var_name}\", it is not in the list of hists to include")
                continue

            fill_info = {
                var_name    : var_values,
                "process"   : hist_axis_name,
                "weight"    : event_weights[event_selection_mask],
                "eft_coeff" : eft_coeffs_cut,
            }

            hout[var_name].fill(**fill_info)

        return hout


    def postprocess(self, accumulator):
        return accumulator
