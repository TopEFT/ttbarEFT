import awkward as ak
import numpy as np 

triggers_dict = {
    "2016APV": {
        "SingleMuon": [
            "Mu50", 
            "TkMu50",
        ],
        "DoubleEG":[
            "DoubleEle33_CaloIdL_MW",
        ],
    },
    "2016APV_run276453": {
        "SingleMuon": [
            "Mu50", 
            "TkMu50",
        ],
        "DoubleEG":[
            "DoubleEle33_CaloIdL_GsfTrkIdVL",
        ],
    },
    "2016": {
        "SingleMuon": [
            "Mu50", 
            "TkMu50",
        ],
        "DoubleEG":[
            "DoubleEle33_CaloIdL_GsfTrkIdVL_MW"
        ],
    },
    "2017": {
        "SingleMuon":[
            "Mu50",
            "TkMu100",
            "OldMu100", 
        ],
        "DoubleEG":[
            "DoubleEle33_CaloIdL_MW",
        ],
    },
    "2018": {
        "SingleMuon":[
            "Mu50", 
            "TkMu100",
            "OldMu100",
        ],
        "EGamma":[
            "DoubleEle25_CaloIdL_MW",
        ],
    },
}

exclude_triggers_dict = {
    "2016": {
        "SingleMuon": [],
        "DoubleEG": triggers_dict["2016"]["SingleMuon"],
    },  
    "2017": {
        "SingleMuon": [],
        "DoubleEG": triggers_dict["2017"]["SingleMuon"],
    },
    "2018": {
        "SingleMuon": [],
        "EGamma": triggers_dict["2018"]["SingleMuon"],
    },
}

MC_lepcat_triggers_dict = {
    "ee": "DoubleEG",   #for 2018, ee trigger is from EGamma
    "em": "SingleMuon", 
    "mm": "SingleMuon",
}

# This is a helper function called by trg_pass_no_overlap
#   - Takes events objects, and a lits of triggers
#   - Returns an array the same length as events, elements are true if the event passed at least one of the triggers and false otherwise
def passes_trg_inlst(events,trg_name_lst):
    tpass = np.zeros_like(np.array(events.MET.pt), dtype=bool)
    trg_info_dict = events.HLT

    # "fields" should be list of all triggers in the dataset
    common_triggers = set(trg_info_dict.fields) & set(trg_name_lst)

    # Check to make sure that at least one of our specified triggers is present in the dataset
    if len(common_triggers) == 0 and len(trg_name_lst):
        raise Exception("No triggers from the sample matched to the ones used in the analysis.")

    for trg_name in common_triggers:
        tpass = tpass | trg_info_dict[trg_name]
    return tpass


# This is what we call from the processor
#   - Returns an array the len of events
#   - Elements are false if they do not pass any of the triggers defined in dataset_dict
#   - In the case of data, events are also false if they overlap with another dataset
def trg_pass_no_overlap(events,is_data,dataset,year,dataset_dict,exclude_dict,lep_cat):

    # Initialize ararys and lists, get trg pass info from events
    trg_passes    = np.zeros_like(np.array(events.MET.pt), dtype=bool) # Array of False the len of events
    trg_overlaps  = np.zeros_like(np.array(events.MET.pt), dtype=bool) # Array of False the len of events
    trg_info_dict = events.HLT
    full_trg_lst  = []

    run_number = events.run

    # In case of data, check if events overlap with other datasets
    if is_data:
        # use different triggers for 2016APV for runs >= 276453
        if (year == '2016APV') and (run_number >= 276453): 
            dataset= "2016APV_run276453"

        trg_passes = passes_trg_inlst(events,dataset_dict[year][dataset])
        trg_overlaps = passes_trg_inlst(events, exclude_dict[year][dataset])

    # In case of MC, pick the trigger dictionary dataset based on lepton channel
    # ee channel is filled from electron triggers, em and mm channels from only muon triggers
    else: 
        if (year == '2018') and (lep_cat == 'ee'):
            dataset_name = 'EGamma'
        else: 
            dataset_name = MC_lepcat_triggers_dict[lep_cat]

        print(f"\n\n lep_cat: {lep_cat} \n dataset_name:{dataset_name} \n dictionary: {dataset_dict[year][dataset_name]} \n\n")
        trg_passes = passes_trg_inlst(events, dataset_dict[year][dataset_name])


    # Return true if passes trg and does not overlap
    return (trg_passes & ~trg_overlaps)


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
    events['is_e'] = ((n_e_2l==1) & (n_m_2l==0))
    events['is_m'] = ((n_e_2l==0) & (n_m_2l==1))

    # 2l masks
    events['is_ee'] = ((n_e_2l==2) & (n_m_2l==0))
    events['is_em'] = ((n_e_2l==1) & (n_m_2l==1))
    events['is_mm'] = ((n_e_2l==0) & (n_m_2l==2))


def addLepSFs(events, ele, mu):
    # TODO: probably need to change these to SF_2l_ee, SF_2l_mm, and add SF_2l_em 
    # where SF_2l_em = leps[0].SF_ele * leps[0].SF_muon * leps[1].SF_ele * leps[1].SF_muon
    
    leps = events.leps_pt_sorted
    padded_leps = ak.pad_none(leps, 2)

    # leps = ak.concatenate([ele, mu], axis=1)
    # padded_leps = ak.pad_none(leps[ak.argsort(leps.pt, axis=-1,ascending=False)], 2) 

    print(f"\n\n {padded_leps.fields} \n\n")

    events['SF_2l_ee'] = padded_leps[:,0].SF_ele_nom * padded_leps[:,1].SF_ele_nom
    events['SF_2l_ee_up'] = padded_leps[:,0].SF_ele_up * padded_leps[:,1].SF_ele_up
    events['SF_2l_ee_down'] = padded_leps[:,0].SF_ele_down * padded_leps[:,1].SF_ele_down

    events['SF_2l_mm'] = padded_leps[:,0].SF_muon_nom * padded_leps[:,1].SF_muon_nom
    events['SF_2l_mm_up'] = padded_leps[:,0].SF_muon_up * padded_leps[:,1].SF_muon_up
    events['SF_2l_mm_down'] = padded_leps[:,0].SF_muon_down * padded_leps[:,1].SF_muon_down

    events['SF_2l_em'] = padded_leps[:,0].SF_ele_nom * padded_leps[:,0].SF_muon_nom * padded_leps[:,1].SF_ele_nom * padded_leps[:,1].SF_muon_nom
    events['SF_2l_em_up'] = padded_leps[:,0].SF_ele_up * padded_leps[:,0].SF_muon_up * padded_leps[:,1].SF_ele_up * padded_leps[:,1].SF_muon_up
    events['SF_2l_em_down'] = padded_leps[:,0].SF_ele_down * padded_leps[:,0].SF_muon_down * padded_leps[:,1].SF_ele_down * padded_leps[:,1].SF_muon_down
