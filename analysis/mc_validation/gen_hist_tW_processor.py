#!/usanr/bin/env python
import awkward as ak
import numpy as np

import hist
from hist import Hist

from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods import vector

from mt2 import mt2
from mt2 import mt2_arxiv

from topcoffea.modules import utils
import topcoffea.modules.eft_helper as efth

np.seterr(divide='ignore', invalid='ignore', over='ignore')
NanoAODSchema.warn_missing_crossrefs = False

SM_pt = {"ctGIm": 0.0, "ctGRe":0.0, "cHQ3": 0.0, "ctWRe": 0.0, "cbWRe": 0.0, "cHtbRe": 0.0,
         "cQl311": 0.0, "cQl322": 0.0, "cQl333": 0.0,
         "cleQt3Re11": 0.0, "cleQt3Re22": 0.0, "cleQt3Re33": 0.0,
         "cleQt1Re11": 0.0, "cleQt1Re22": 0.0, "cleQt1Re33": 0.0}

# starting point of tWLO_rwgt1 sample
rwgt1 = {"ctGIm": -0.2, "ctGRe":-0.2, "cHQ3": 1.5, "ctWRe": -1.0, "cbWRe": -10.0, "cHtbRe": 4.0,
         "cQl311": 4.0, "cQl322": 4.1, "cQl333": 4.2,
         "cleQt3Re11": 12.0, "cleQt3Re22": 12.1, "cleQt3Re33": 12.2,
         "cleQt1Re11": 15.0, "cleQt1Re22": 15.1, "cleQt1Re33": 15.2}

# starting point of tWLO_rwgt2 sample
rwgt2 = {"ctGIm": -0.5, "ctGRe":-0.5, "cHQ3": 1.5, "ctWRe": -1.0, "cbWRe": -8.0, "cHtbRe": 4.0,
         "cQl311": 8.0, "cQl322": 8.0, "cQl333": 8.0,
         "cleQt3Re11": 12.0, "cleQt3Re22": 12.0, "cleQt3Re33": 12.0,
         "cleQt1Re11": 15.0, "cleQt1Re22": 15.0, "cleQt1Re33": 15.0}

rwgt_choice = rwgt2

# Get the lumi for the given year (in fb^-1)
# Ref: https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2#Luminosity_for_pp_13_TeV_data
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

# Create list of wc values in the correct order
def order_wc_values(wcs, ref_pts):
    '''Returns list of wc values in the same order as the list of wc names based on a dictionary
    '''
    wc_names = wcs
    ref_pts = ref_pts

    wc_values = []
    for i in wc_names:
        wc_values.append(ref_pts[i])

    return wc_values

# Calculate event weights from wc values and eft fit coefficients
def calc_event_weights(eft_coeffs, wc_vals):
    '''Returns an array that contains the event weight for each event.
    eft_coeffs: Array of eft fit coefficients for each event
    wc_vals: wilson coefficient values desired for the event weight calculation, listed in the same order as the wc_lst
             such that the multiplication with eft_coeffs is correct
             The correct ordering can be achieved with the order_wc_values function
    '''
    event_weight = np.empty_like(eft_coeffs)

    wcs = np.hstack((np.ones(1),wc_vals))
    wc_cross_terms = []
    index = 0
    for j in range(len(wcs)):
        for k in range (j+1):
            term = wcs[j]*wcs[k]
            wc_cross_terms.append(term)
    event_weight = np.sum(np.multiply(wc_cross_terms, eft_coeffs), axis=1)

    return event_weight

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

    # mt2_val = mt2(...)
    mt2_val = mt2(
        l0.mass, l0.px, l0.py,                          # visible particle #1
        l1.mass, l1.px, l1.py,                          # visible particle #2 
        misspart.px, misspart.py,                       # missing transverse momentum
        np.zeros_like(met.pt), np.zeros_like(met.pt)    # invisible masses
    )

    return mt2_val

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst = None, dtype=np.float32, do_errors=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst

        self._dtype = dtype
        self._do_errors = do_errors

        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)

        # Create the histograms with new scikit hist
        self._histo_dict = {
            "njets"             : Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="njets", label='njets'), storage="weight"),
            "nleps"             : Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="nleps", label='nleps'), storage="weight"),
            "ntops"             : Hist(hist.axis.Regular(bins=10, start=0, stop=10, name="ntops", label='ntops'), storage="weight"),
            "mll"               : Hist(hist.axis.Regular(bins=40, start=0, stop=800, name='mll', label='invariant mass of leptons [GeV]'), storage="weight"),
            "dr_leps"           : Hist(hist.axis.Regular(bins=24, start=0, stop=6, name='dr_leps', label="$\Delta R$ (leading lepton, subleading lepton)"), storage="weight"),
            "l0pt"         		: Hist(hist.axis.Regular(bins=40, start=0, stop=400, name='l0pt', label="leading lepton $p_T$ [GeV]"), storage="weight"),
            "top_pt"            : Hist(hist.axis.Regular(bins=70, start=0, stop=700, name="top_pt", label="$p_T$ of the top quark [GeV]"), storage="weight"),
            "mt2"               : Hist(hist.axis.Variable([0,10,50,80,100,110,120,140,160,180,220,500], name="mt2", label="$m_{T2}$ [GeV]"), storage="weight"),
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

        # Standard 2l, 1j selections
        at_least_two_leps = ak.fill_none(nleps>=2,False)
        at_least_one_jet = ak.fill_none(njets>=1,False)
        at_least_one_top = ak.fill_none(ntops==1, False)

        selections = PackedSelection()
        selections.add('2l', at_least_two_leps)
        selections.add('1j', at_least_one_jet)
        selections.add('1t', at_least_one_top)
        event_selection_mask = selections.all('2l', '1j', '1t')

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


        eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

        if eft_coeffs is None:
            genw = events["genWeight"]
            event_weights = genw[event_selection_mask]
        else: 
            wc_lst = order_wc_values(self._wc_names_lst, rwgt_choice)
            event_weights = calc_event_weights(eft_coeffs_cut, wc_lst)

        ######## Normalization ########

        # Normalize by (xsec/sow)
        lumi = 1000.0*get_lumi(year) #lumi is fb^-1, but xsec is in pb

        # tW_yukawa_const = 1.5217
        # tW_new1_const = 1.5007769740968973
        tW_powheg_const = 1.845
        tWLO_stpt1_const = 1.87421
        tWLO_stpt2_const = 1.88947


        norm = (xsec/sow)*lumi
        # norm = (xsec/sow)*(1/tW_yukawa_const)*lumi
        # norm = (xsec/sow)*(1/tW_new1_const)
        # norm = (xsec/sow)*(1/tW_powheg_const)

        # norm = (xsec/sow)*(1/tWLO_stpt1_const)
        # norm = (xsec/sow)*(1/tWLO_stpt2_const)
        # norm = (xsec/sow)*(1/tW_powheg_const)

        counts = np.ones_like(events['event'])[event_selection_mask]

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

        for var_name, var_values in variables_to_fill.items():
            if var_name not in self._hist_lst:
                print(f"Skipping \"{var_name}\", it is not in the list of hists to include")
                continue

            hout[var_name].fill(var_values, weight=event_weights*norm)

        return hout

    def postprocess(self, accumulator):
        return accumulator

