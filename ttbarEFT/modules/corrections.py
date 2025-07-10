import numpy as np
import pickle
import gzip
import awkward as ak

import correctionlib

from coffea import lookup_tools

from topcoffea.modules.paths import topcoffea_path
from ttbarEFT.modules.paths import ttbarEFT_path

def ApplyJetVetoMaps(jets, year):

    if year == "2016APV": 
        fname = ttbarEFT_path("data/POG/JME/2016preVFP_UL/jetvetomaps.json.gz")
        key = "Summer19UL16_V1"
    elif year == "2016": 
        fname = ttbarEFT_path("data/POG/JME/2016postVFP_UL/jetvetomaps.json.gz")
        key = "Summer19UL16_V1"
    elif year == "2017": 
        fname = ttbarEFT_path("data/POG/JME/2017_UL/jetvetomaps.json.gz")
        key = "Summer19UL17_V1"
    elif year == "2018":
        fname = ttbarEFT_path("data/POG/JME/2018_UL/jetvetomaps.json.gz")
        key = "Summer19UL18_V1"

    # Grab the json
    ceval = correctionlib.CorrectionSet.from_file(fname)

    # Flatten the inputs
    eta_flat = ak.flatten(jets.eta)
    phi_flat = ak.flatten(jets.phi)

    #Put mins and maxes on the accepted values
    eta_flat_bound = ak.where(eta_flat>5.19,5.19,ak.where(eta_flat<-5.19,-5.19,eta_flat))
    phi_flat_bound = ak.where(phi_flat>3.14159,3.14159,ak.where(phi_flat<-3.14159,-3.14159,phi_flat))

    #Get pass/fail values for each jet (0 is pass and >0 is fail)
    jet_vetomap_flat = ceval[key].evaluate('jetvetomap',eta_flat_bound,phi_flat_bound)
    
    #Unflatten the array
    jet_vetomap_score = ak.unflatten(jet_vetomap_flat,ak.num(jets.phi))

    #Sum the outputs for each event (if the sum is >0, the event will fail)
    veto_map_event = ak.sum(jet_vetomap_score, axis=-1)

    return veto_map_event

