import awkward as ak

triggers_dict = {
    "2016": {
        "SingleElectron":[
            "Ele27_WPTight_Gsf",
            "Photon175",    
        ],
        "SingleMuon": [
            "Mu50", 
            "TkMu50",
        ],
        "DoubleEG":[
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ],
    },
    "2017": {
        "SingleElectron":[
            "Ele35_WPTight_Gsf",
            "Photon200",
        ],
        "SingleMuon":[
            "Mu50",
            "TkMu100",
            "OldMu100", 
        ],
        "DoubleEG":[
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        ],
    },
    "2018": {
        "EGamma":[
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Ele32_WPTight_Gsf",
            "Photon200",
        ],
        "SingleMuon":[
            "Mu50", 
            "TkMu100",
            "OldMu100",
        ],
    },
}

exclude_triggers_dict = {
    "2016": {
        "SingleMuon": [],
        "SingleElectron": triggers_dict["2016"]["SingleMuon"],
        "DoubleEG": triggers_dict["2016"]["SingleMuon"] + triggers_dict["2016"]["SingleElectron"],
    },  
    "2017": {
        "SingleMuon": [],
        "SingleElectron": triggers_dict["2017"]["SingleMuon"],
        "DoubleEG": triggers_dict["2017"]["SingleMuon"] + triggers_dict["2017"]["SingleElectron"],
    },
    "2018": {
        "SingleMuon": [],
        "EGamma": triggers_dict["2018"]["SingleMuon"],
    },
}


def addLepCatMasks(events):
    
    leps = events.leps_pt_sorted

    padded_leps = ak.pad_none(leps, 2)
    padded_leps_id = padded_leps.pdgId

    is_e_mask = (abs(padded_leps_id)==11)
    is_m_mask = (abs(padded_leps_id)==13)

    n_e_1l = ak.sum(is_e_mask[:,0:1],axis=-1) # Make sure we only look at first lep
    n_m_1l = ak.sum(is_m_mask[:,0:1],axis=-1) # Make sure we only look at first lep
    n_e_2l = ak.sum(is_e_mask[:,0:2],axis=-1) # Make sure we only look at first two leps
    n_m_2l = ak.sum(is_m_mask[:,0:2],axis=-1) # Make sure we only look at first two leps        


    # 1l masks
    events["is_e"] = ((n_e_2l==1) & (n_m_2l==0))
    events["is_m"] = ((n_e_2l==0) & (n_m_2l==1))

    # 2l masks
    events['is_ee'] = ((n_e_2l==2) & (n_m_2l==0))
    events['is_em'] = ((n_e_2l==1) & (n_m_2l==1))
    events['is_mm'] = ((n_e_2l==0) & (n_m_2l==2))

def addLepSFs(events):

    leps = events.leps_pt_sorted
    padded_leps = ak.pad_none(leps, 2)

    print(f"\n\n padded_leps.fields: {padded_leps.fields} \n\n")

    events['SF_2l_ele'] = padded_leps[:,0].SF_ele_nom * padded_leps[:,1].SF_ele_nom
    events['SF_2l_ele_up'] = padded_leps[:,0].SF_ele_up * padded_leps[:,1].SF_ele_up
    events['SF_2l_ele_down'] = padded_leps[:,0].SF_ele_down * padded_leps[:,1].SF_ele_down

    events['SF_2l_muon'] = padded_leps[:,0].SF_muon_nom * padded_leps[:,1].SF_muon_nom
    events['SF_2l_muon_up'] = padded_leps[:,0].SF_muon_up * padded_leps[:,1].SF_muon_up
    events['SF_2l_muon_down'] = padded_leps[:,0].SF_muon_down * padded_leps[:,1].SF_muon_down

