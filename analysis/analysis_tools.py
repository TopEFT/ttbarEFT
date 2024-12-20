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
    nu    = leps[(abs(leps.pdgId) == 12) | (abs(leps.pdgId) == 14)]

    ######## Lep selection ########
    el    = el[(el.pt>20) & (abs(el.eta)<2.5)]
    mu    = mu[(mu.pt>20) & (abs(mu.eta)<2.5)]
    leps  = ak.concatenate([el,mu],axis=1)

    ######## Jet selection ########
    jets  = jets[(jets.pt>30) & (abs(jets.eta)<2.5)]
    jets  = jets[isClean(jets, leps, drmin=0.4) & isClean(jets, nu, drmin=0.4)]

    return leps, jets

def genEventSelection(leps, jets):
    nleps = ak.num(leps)
    njets = ak.num(jets)

    at_least_two_leps = ak.fill_none(nleps>=2,False)
    at_least_two_jets = ak.fill_none(njets>=2, False)
    
    selections = PackedSelection()
    selections.add('2l', at_least_two_leps)
    selections.add('2j', at_least_two_jets)
    event_selection_mask = selections.all('2l', '2j')
    
    return event_selection_mask