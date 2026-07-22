from coffea.analysis_tools import PackedSelection
from coffea.processor import AccumulatorABC

import awkward as ak
import numpy as np
import torch

class TensorAccumulator(AccumulatorABC):
    def __init__(self, tensor: torch.Tensor, dtype=torch.float64):
        self._tensor = tensor
        self._dtype = dtype
        
    def add(self, other: "TensorAccumulator") -> "TensorAccumulator":
        return TensorAccumulator(torch.concat([self._tensor, other._tensor]))
    
    def __add__(self, other):
        return self.add(other)

    def __iadd__(self, other):
        return self.add(other)
        
    def get(self) -> torch.Tensor:
        return self._tensor

    def identity(self):
        return TensorAccumulator(torch.Tensor)
        
    def concat(self, tensor: torch.Tensor):
        return TensorAccumulator(torch.concat([self._tensor, tensor], axis=0).type(self._dtype))


def isClean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)

def genObjectSelection(events):
    ######## Initialize objets ########
    leps  = events.GenDressedLepton
    jets  = events.GenJet
    el    = leps[abs(leps.pdgId) == 11]
    mu    = leps[abs(leps.pdgId) == 13]

    ######## Lep selection ########
    el    = el[(el.pt>20) & (abs(el.eta)<2.5)]
    mu    = mu[(mu.pt>20) & (abs(mu.eta)<2.5)]
    leps  = ak.concatenate([el,mu],axis=1)
    leps  = leps[ak.argsort(leps.pt, axis=-1, ascending=False)]

    ######## Jet selection ########
    jets  = jets[(jets.pt>30) & (abs(jets.eta)<2.5)]
    jets  = jets[isClean(jets, leps, drmin=0.4)]
    jets  = jets[ak.argsort(jets.pt, axis=-1, ascending=False)]

    return leps, jets

def genEventSelection(leps, jets):
    pMask = leps.pdgId > 0
    nMask = leps.pdgId < 0
    njets = ak.num(jets)
    
    os_leps = ak.fill_none(ak.any(pMask, 1) & ak.any(nMask, 1), False)
    at_least_two_jets = ak.fill_none(njets>=2, False)

    selections = PackedSelection()
    selections.add('osl', os_leps)
    selections.add('2j', at_least_two_jets)
    event_selection_mask = selections.all('osl', '2j')
    
    return event_selection_mask

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

    # switch to mt2_val = mt2(...) after updating mt2 version
    mt2_val = mt2_arxiv(
        l0.mass, l0.px, l0.py,                          # visible particle #1
        l1.mass, l1.px, l1.py,                          # visible particle #2 
        misspart.px, misspart.py,                       # missing transverse momentum
        np.zeros_like(met.pt), np.zeros_like(met.pt)    # invisible masses
    )

    return mt2_val
