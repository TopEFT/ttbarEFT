import numpy as np
import awkward as ak

class Run2LeptonSelection:
    def __init__(self):
        pass

    def is_sel_ele(self, ele):
        # apply eta, HEEP cuts
        pt_mask = ele.pt > 40 # orig 35, inc to 40 b/c of trigger eff validity
        eta_mask = ((abs(ele.eta) < 1.4442) | (abs(ele.eta) > 1.556)) & (abs(ele.eta) < 2.4)
        HEEP_mask = ele.cutBased_HEEP #bool, a simple robust ID designed to be safe for high electrons

        return (pt_mask & eta_mask & HEEP_mask)

    def is_sel_muon(self, muon):
        # apply eta, highPtId and tkIsoId cuts
        pt_mask = muon.pt > 53
        eta_mask = (abs(muon.eta) < 2.4)
        highPtId_mask = (muon.highPtId == 2) #global high pT, includes tracker high pT
        tkIsoId_mask = (muon.tkIsoId > 0) #selects tkIsoLoose & tkIsoTight

        return (pt_mask & eta_mask & highPtId_mask & tkIsoId_mask)


def is_pres_jet(jet):

    pt_mask = jet.pt > 30
    eta_mask = (abs(jet.eta) < 2.4)
    id_mask = jet.jetId >= 6 #2016: jetID==7 means pass loose, tight, tightLepVeto ID; 2017-2018: jetID==6 means pass tight and tightLepVeto ID
    # puId_mask = jet.puId > 1 
    puId_mask = (jet.puId > 1) | (jet.pt >= 50) # Jets with pT>=50 always pass 

    return (pt_mask & eta_mask & id_mask & puId_mask)


def isClean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)