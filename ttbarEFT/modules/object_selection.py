import numpy as np
import awkward as ak

class Run2LeptonSelection:
    def __init__(self):
        pass

    def is_sel_ele(self, ele):
        # apply eta, HEEP cuts
        pt_mask = ele.pt > 35
        eta_mask = ((abs(ele.eta) < 1.4442) | (abs(ele.eta) > 1.556)) & (abs(ele.eta) < 2.4)
        HEEP_mask = ele.cutBased_HEEP

        return (pt_mask & eta_mask & HEEP_mask)

    def is_sel_muon(self, muon):
        # apply eta, highPtId and tkIsoId cuts
        pt_mask = muon.pt > 53
        eta_mask = (abs(muon.eta) < 2.4)
        highPtId_mask = (muon.highPtId == 2)
        tkIsoId_mask = (muon.tkIsoId > 0)

        return (pt_mask & eta_mask & highPtId_mask & tkIsoId_mask)


def is_pres_jet(jet):

    eta_mask = (abs(jet.eta) < 2.4)
    id_mask = jet.jetId >= 6 #using jet.isTightLeptonVeto==True should be equivalent
    puId_mask = jet.puId > 1 
    pt_mask = jet.pt > 30

    return (id_mask & puId_mask & pt_mask & eta_mask)


def isClean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)