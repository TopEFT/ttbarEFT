import numpy as np
import awkward as ak

def isClean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)

def gen_lepton_selections(ele, mu):
	e_selec = ((ele.pt>20) & (abs(ele.eta)<2.5))
	m_selec = ((mu.pt>20) & (abs(mu.eta)<2.5))

	leps = ak.concatenate([ele[e_selec],mu[m_selec]],axis=1)	#apply selections
	leps = leps[ak.argsort(leps.pt, axis=-1, ascending=False)]	#sort leps by pt

	return leps

def gen_jets_selection(jets, leps, nu):
	j_selec = ((jets.pt>30) & (abs(jets.eta)<2.5))
	jets = jets[j_selec]										#apply selections

	jets= jets[isClean(jets, leps, drmin=0.4) & isClean(jets, nu, drmin=0.4)]	# clean the jets
	jets = jets[ak.argsort(jets.pt, axis=-1, ascending=False)]	# sort jets in pt

	return jets 
	