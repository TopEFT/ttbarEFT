import awkward as ak
import numpy as np 
import copy

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
    "2016APV_run276453": {
        "SingleMuon": [],
        "DoubleEG": triggers_dict["2016APV_run276453"]["SingleMuon"],
    },  
    "2016APV": {
        "SingleMuon": [],
        "DoubleEG": triggers_dict["2016APV"]["SingleMuon"],
    },  
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
    "ee": "DoubleEG",   #for 2018, ee trigger is from EGamma, hardcoded in trigger function
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
        # if (year == '2016APV') and (run_number >= 276453): 
        #     dataset= "2016APV_run276453"
        # trg_passes = passes_trg_inlst(events,dataset_dict[year][dataset])
        # trg_overlaps = passes_trg_inlst(events, exclude_dict[year][dataset])

        if (year == '2016APV'): 
            run_mask = (run_number >= 276453)

            # Process Early
            if ak.any(run_mask):
                trg_passes = ak.where(run_mask, passes_trg_inlst(events, dataset_dict["2016APV_run276453"][dataset]), trg_passes)
                trg_overlaps = ak.where(run_mask, passes_trg_inlst(events, exclude_dict["2016APV_run276453"][dataset]), trg_overlaps)

            # Process Late
            if ak.any(~run_mask):
                trg_passes = ak.where(~run_mask, passes_trg_inlst(events, dataset_dict["2016APV"][dataset]), trg_passes)
                trg_overlaps = ak.where(~run_mask, passes_trg_inlst(events, exclude_dict["2016APV"][dataset]), trg_overlaps)

        else: 
            trg_passes = passes_trg_inlst(events,dataset_dict[year][dataset])
            trg_overlaps = passes_trg_inlst(events, exclude_dict[year][dataset])

    # In case of MC, pick the trigger dictionary dataset based on lepton channel
    # ee channel is filled from electron triggers, em and mm channels from only muon triggers
    else: 
        if (year == '2018') and (lep_cat == 'ee'):
            dataset_name = 'EGamma'
        else: 
            dataset_name = MC_lepcat_triggers_dict[lep_cat]

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

def add2losMask(events, year, isData):
    
    leps = events.leps_pt_sorted
    padded_leps = ak.pad_none(leps,2)

    filter_flags = events.Flag
    filters = filter_flags.goodVertices & filter_flags.globalSuperTightHalo2016Filter & filter_flags.HBHENoiseFilter & filter_flags.HBHENoiseIsoFilter & filter_flags.EcalDeadCellTriggerPrimitiveFilter & filter_flags.BadPFMuonFilter & filter_flags.BadPFMuonDzFilter & (((year == "2016")|(year == "2016APV")) | filter_flags.ecalBadCalibFilter) & filter_flags.eeBadScFilter # (isData | filter_flags.eeBadScFilter)

    llpairs = ak.combinations(padded_leps, 2, fields=["l0","l1"])
    events["minMllAFAS"] = ak.min((llpairs.l0+llpairs.l1).mass, axis=-1)
    # cleanup = events.minMllAFAS > 12
    cleanup = events.minMllAFAS > 20

    dilep = (ak.num(leps)) >= 2
    exclusive = (ak.num(leps)) < 3

    os = (padded_leps[:,0].charge + padded_leps[:,1].charge)==0
    mask = (filters & cleanup & dilep & exclusive & os)

    events['is2los'] = ak.fill_none(mask,False)

def addmllMasks(events):

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]

    mll = (l0+l1).mass
    
    mll_mask = mll > 106
    events['mllg106'] = ak.fill_none(mll_mask, False)

    DY_mass_mask = ((mll < 75) | (mll > 105))
    events['mllDYmask'] = ak.fill_none(DY_mass_mask, False)


def getHemMask(events, year, isData, jets): 
    '''
    Returns: 
        ~is_hem_events (bool mask), used for coffea PackedSelection. True for events that should pass. False for events with HEM lep/jet
        hem_wieght (float), event weight = 0.35 for events with HEM lep/jet. else = 1.0 (no modification)
    '''

    hem_mask = ak.zeros_like(events.event, dtype=bool)      # True if contains HEM event 
    hem_weight = ak.ones_like(events.event, dtype=float)        # event weight multiplier

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]

    if year == '2018': 
        hem_jets = (jets.eta < -1.3) & (jets.phi < -1.57) & (jets.phi > -0.87)
        event_has_hem_jet = ak.any(hem_jets, axis=1)


        l0_is_hem_ele = (abs(l0.pdgId) == 11) & (l0.eta < -1.3) & (l0.phi < -0.87) & (l0.phi > -1.57)
        l1_is_hem_ele = (abs(l1.pdgId) == 11) & (l1.eta < -1.3) & (l1.phi < -0.87) & (l1.phi > -1.57)
        event_has_hem_lep = ak.fill_none(l0_is_hem_ele | l1_is_hem_ele, False)

        is_hem_event = event_has_hem_jet | event_has_hem_lep

        if isData: 
            hem_mask = is_hem_event & (events.run >= 319077)
        else: 
            hem_weight = ak.where(is_hem_event, 0.35, 1.0)          # if it's a hem event, give a weight of 0.35
            # is_hem_event = ak.fill_none(ak.zeros_like(is_hem_event, dtype=bool), False)

    return ~hem_mask, hem_weight
