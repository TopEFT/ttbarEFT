import numpy as np
import pickle
import gzip
import awkward as ak

import correctionlib

from coffea import lookup_tools

from topcoffea.modules.paths import topcoffea_path
from ttbarEFT.modules.paths import ttbarEFT_path

clib_year_map = {
    "2016APV": "2016preVFP_UL",
    "2016preVFP": "2016preVFP_UL",
    "2016": "2016postVFP_UL",
    "2017": "2017_UL",
    "2018": "2018_UL",
    "2022": "2022_Summer22",
    "2022EE": "2022_Summer22EE",
    "2023": "2023_Summer23",
    "2023BPix": "2023_Summer23BPix",
}

def ApplyJetVetoMaps(jets, year):

    jet_veto_dict = {
        "2016APV": "Summer19UL16_V1",
        "2016": "Summer19UL16_V1",
        "2017": "Summer19UL17_V1",
        "2018": "Summer19UL18_V1",
        "2022": "Summer22_23Sep2023_RunCD_V1",
        "2022EE": "Summer22EE_23Sep2023_RunEFG_V1",
        "2023": "Summer23Prompt23_RunC_V1",
        "2023BPix": "Summer23BPixPrompt23_RunD_V1"
    }

    jme_year = clib_year_map[year]
    key = jet_veto_dict[year]
    json_path = ttbarEFT_path(f"data/POG/JME/{jme_year}/jetvetomaps.json.gz")

    # Grab the json
    ceval = correctionlib.CorrectionSet.from_file(json_path)

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
