import numpy as np
import pickle
import gzip
import awkward as ak

import correctionlib

from coffea import lookup_tools
from coffea.lookup_tools import txt_converters, rochester_lookup

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

def AttachElectronTrigSF(electrons, year):

    extLepSF = lookup_tools.extractor()

    extLepSF.add_weight_sets([f"ElecSF_2016APV_barrel UL2016preVFP_Barrel_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])
    extLepSF.add_weight_sets([f"ElecSF_2016APV_endcap UL2016preVFP_Endcaps_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])
    extLepSF.add_weight_sets([f"ElecSF_2016_barrel UL2016postVFP_Barrel_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])
    extLepSF.add_weight_sets([f"ElecSF_2016_endcap UL2016postVFP_Endcaps_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])
    extLepSF.add_weight_sets([f"ElecSF_2017_barrel UL2017_Barrel_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])
    extLepSF.add_weight_sets([f"ElecSF_2017_endcap UL2017_Endcaps_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])
    extLepSF.add_weight_sets([f"ElecSF_2018_barrel UL2018_Barrel_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])
    extLepSF.add_weight_sets([f"ElecSF_2018_endcap UL2018_Endcaps_Et {ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')}"])

    extLepSF.finalize()
    SFevaluator = extLepSF.make_evaluator()

    eta = electrons.eta
    eta_barrel_mask = ak.flatten(abs(eta) < 1.4442)
    
    pt = electrons.pt

    ElecSF = ak.where(
                eta_barrel_mask, 
                SFevaluator[f"ElecSF_{year}_barrel"](pt, eta),
                SFevaluator[f"ElecSF_{year}_endcap"](pt, eta)
            )

    electrons['SF_elec_trig_nom'] = ElecSF 


def AttachElectronSF(electrons, year):
    '''
    Inserts the following scale factor values inot the electrons array 
        'sf_nom': nominal 
        'sf_up': up
        'sf_down': down

    In Run2, the pT regions are above and below 20 GeV. We only selection ele pT>20 so we only 
    have to apply the RecoAbove20 region 
    To include Run3, this function would need to be updated for Run2 vs Run3 since there are 3 pT regions
    See topeft/topeft/modules/corrections.py AttachElectronSF as an example 
    '''

    if year not in clib_year_map.keys():
        raise Exception(f"Error: Unknown year \"{year}\".")

    egm_tag_map = {
        "2016preVFP_UL": "2016preVFP",
        "2016postVFP_UL": "2016postVFP",
        "2017_UL": "2017",
        "2018_UL": "2018",
        "2022_Summer22": "2022Re-recoBCD",
        "2022_Summer22EE": "2022Re-recoE+PromptFG",
        "2022_Summer23": "2023PromptC",
        "2022_Summer23BPix": "2023PromptD",
    }

    HEEPSF_B = {
        "2016preVFP" : 0.985,
        "2016postVFP" : 0.985,
        "2017" : 0.979,
        "2018" : 0.973,
    }

    HEEPSF_E = {
        "2016preVFP" : 0.990,
        "2016postVFP" : 0.990,
        "2017" : 0.987,
        "2018" : 0.980,   
    }

    # initialize electron variables
    eta = electrons.eta
    pt = electrons.pt
    phi = electrons.phi

    eta_flat = ak.flatten(eta)
    pt_flat = ak.flatten(pt)
    phi_flat = ak.flatten(phi)

    # initilize clib electron SFs file
    clib_year = clib_year_map[year]
    egm_year = egm_tag_map[clib_year]
    egm_tag = "Electron-ID-SF"
    egm_tag = "UL-" + "Electron-ID-SF"
    # json_path = topcoffea_path(f"data/POG/EGM/{clib_year}/electron.json.gz")
    json_path = ttbarEFT_path(f"data/POG/EGM/{clib_year}/electron.json.gz")
    ceval = correctionlib.CorrectionSet.from_file(json_path)

    # create the pt mask
    pt_tag = "RecoAbove20" 
    pt_mask = ak.flatten(pt >= 20)
    pt_mod = ak.where(pt_mask, pt_flat, 900) #returns (pt_val=pt_flat if pt_mask=True), else pt_val=900


    ###### evaluate reco sf ######

    egm_args = [pt_tag, eta_flat, pt_mod]
    # arrays where if pt_mask=True, get sf from clib; else = 1
    reco_nom_flat = ak.to_numpy(ak.where(pt_mask, ceval[egm_tag].evaluate(egm_year, "sf", *egm_args), 1))
    reco_up_flat = ak.to_numpy(ak.where(pt_mask, ceval[egm_tag].evaluate(egm_year, "sfup", *egm_args), 1))
    reco_down_flat = ak.to_numpy(ak.where(pt_mask, ceval[egm_tag].evaluate(egm_year, "sfdown", *egm_args), 1))

    reco_nom = ak.unflatten(reco_nom_flat, ak.num(pt))
    reco_up = ak.unflatten(reco_up_flat, ak.num(pt))
    reco_down = ak.unflatten(reco_down_flat, ak.num(pt))


    ###### evaluate HEEP sf ######

    eta_barrel_mask = ak.flatten(abs(eta) < 1.4442)

    HEEP_eta_bins = {
        'barrel': eta_barrel_mask,
        'endcap': ~eta_barrel_mask
    }

    HEEP_pt_bins = {
        'barrel': {
            'less100': [0, 100], 
            'greater100': [100, np.inf]        
        },
        'endcap': {
            'less100': [0, 100],
            '100to300': [100, 300],
            'greater300': [300, np.inf]
        } 
    }

    HEEP_pt_bins_scaling = {
        'barrel': {
            'less100': {
                'up':1.01, 
                'down': 0.99},
            'greater100': {
                'up': (1 + (0.01 + (0.02/900)*(pt_flat-100))),
                'down': (1 - (0.01 + (0.02/900)*(pt_flat-100)))}
        },
        'endcap': {
            'less100': {
                'up': 1.02,
                'down': 0.98
            },
            '100to300': {
                'up': (1 + (0.02 + (0.03/200)*(pt_flat-100))),
                'down': (1 - (0.02 + (0.03/200)*(pt_flat-100))) 
            },
            'greater300': {
                'up': 1.05,
                'down':0.95 
            },
        } 
    }

    # nominal HEEP SF are just based on eta mask for EB/EE
    HEEPSF_nom_flat = ak.where(
            eta_barrel_mask,
            HEEPSF_B[egm_year],
            HEEPSF_E[egm_year]
        )

    # up/down HEEP SF are also binned in pt 
    HEEPSF_up_perbin = []
    HEEPSF_down_perbin = []

    for eta_region, eta_mask in HEEP_eta_bins.items():
        for pt_region, bin_edges in HEEP_pt_bins[eta_region].items(): 
            # eta mask also needs to be applied
            pt_mask = ak.flatten((pt >= bin_edges[0]) & (pt < bin_edges[1]))

            HEEPSF_up_perbin.append(ak.where(
                eta_mask & pt_mask, 
                HEEP_pt_bins_scaling[eta_region][pt_region]['up'], 
                1))

            HEEPSF_down_perbin.append(ak.where(
                eta_mask & pt_mask, 
                HEEP_pt_bins_scaling[eta_region][pt_region]['down'],
                1)) 


    HEEPSF_up_flat = []
    HEEPSF_down_flat = []

    for idr, HEEPSF_up_bin_flat in enumerate(HEEPSF_up_perbin):
        HEEPSF_up_bin_flat = ak.to_numpy(HEEPSF_up_bin_flat)
        HEEPSF_down_bin_flat = ak.to_numpy(HEEPSF_down_perbin[idr])

        if idr==0: 
            HEEPSF_up_flat = HEEPSF_up_bin_flat
            HEEPSF_down_flat = HEEPSF_down_bin_flat
        else: 
            HEEPSF_up_flat = HEEPSF_up_flat * HEEPSF_up_bin_flat
            HEEPSF_down_flat = HEEPSF_down_flat * HEEPSF_down_bin_flat

    HEEPSF_nom = ak.unflatten(HEEPSF_nom_flat, ak.num(eta))
    HEEPSF_up = ak.unflatten(HEEPSF_up_flat, ak.num(eta))
    HEEPSF_down = ak.unflatten(HEEPSF_down_flat, ak.num(eta))

    # Attach SFs (reco*HEEP) to electrons
    electrons['SF_ele_nom'] = reco_nom * HEEPSF_nom
    electrons['SF_ele_up'] = reco_up * HEEPSF_up
    electrons['SF_ele_down'] = reco_up * HEEPSF_down

    electrons['SF_muon_nom'] = ak.ones_like(reco_nom)
    electrons['SF_muon_up'] = ak.ones_like(reco_nom)
    electrons['SF_muon_down'] = ak.ones_like(reco_nom)


def AttachMuonTrigSF(muons, year):
    if year not in clib_year_map.keys():
        raise Exception(f"Error: Unknown year \"{year}\".")

    # initialize muon variables 
    abs_eta = np.abs(muons.eta) #For run3 abs(eta) should be changed to signed eta
    pt = muons.pt 
    phi = muons.phi 

    abseta_flat = ak.flatten(abs_eta)
    pt_flat = ak.flatten(pt)
    phi_flat = ak.flatten(phi) 

    clib_year = clib_year_map[year]
    json_path_HighPt = ttbarEFT_path(f"data/POG/MUO/{clib_year}/muon_HighPt.json.gz")
    json_path_z = ttbarEFT_path(f"data/POG/MUO/{clib_year}/muon_Z.json.gz")

    ceval_HighPt = correctionlib.CorrectionSet.from_file(json_path_HighPt)
    ceval_z = correctionlib.CorrectionSet.from_file(json_path_z)

    trigger_norm_flat = ceval_HighPt["NUM_HLT_DEN_HighPtLooseRelIsoProbes"].evaluate(abseta_flat, pt_flat, "nominal")
    trigger_up_flat = ceval_HighPt["NUM_HLT_DEN_HighPtLooseRelIsoProbes"].evaluate(abseta_flat, pt_flat, "systup")
    trigger_down_flat = ceval_HighPt["NUM_HLT_DEN_HighPtLooseRelIsoProbes"].evaluate(abseta_flat, pt_flat, "systdown")

    trigger_norm = ak.unflatten(trigger_norm_flat, ak.num(pt))
    trigger_up = ak.unflatten(trigger_up_flat, ak.num(pt))
    trigger_down = ak.unflatten(trigger_down_flat, ak.num(pt))

    muons['SF_muon_trig_nom'] = trigger_norm
    muons['SF_muon_trig_up'] = trigger_up
    muons['SF_muon_trig_down'] = trigger_down


def AttachMuonSF(muons, year): 
    '''
    Inserts the following scale factor values into the electrons array 
        'sf_nom': nominal 
        'sf_up': up
        'sf_down': down

    In Run2, the pT regions are above and below 20 GeV. We only selection ele pT>20 so we only 
    have to apply the RecoAbove20 region 
    To include Run3, this function would need to be updated for Run2 vs Run3 since there are 3 pT regions
    See topeft/topeft/modules/corrections.py AttachElectronSF as an example 
    '''

    if year not in clib_year_map.keys():
        raise Exception(f"Error: Unknown year \"{year}\".")

    # initialize muon variables 
    abs_eta = np.abs(muons.eta) #For run3 abs(eta) should be changed to signed eta
    pt = muons.pt 
    phi = muons.phi 

    abseta_flat = ak.flatten(abs_eta)
    pt_flat = ak.flatten(pt)
    phi_flat = ak.flatten(phi)

    clib_year = clib_year_map[year]
    json_path_HighPt = ttbarEFT_path(f"data/POG/MUO/{clib_year}/muon_HighPt.json.gz")
    json_path_z = ttbarEFT_path(f"data/POG/MUO/{clib_year}/muon_Z.json.gz")

    ceval_HighPt = correctionlib.CorrectionSet.from_file(json_path_HighPt)
    ceval_z = correctionlib.CorrectionSet.from_file(json_path_z)

    # pt_mask = ak.flatten(pt >= 50) # the lowest pT bin in muon_HighPt is 50 
    # evaluate SFs

    tracking_nom_flat = ceval_z["NUM_TrackerMuons_DEN_genTracks"].evaluate(abseta_flat, pt_flat, "nominal")
    tracking_up_flat = ceval_z["NUM_TrackerMuons_DEN_genTracks"].evaluate(abseta_flat, pt_flat, "systup")
    tracking_down_flat = ceval_z["NUM_TrackerMuons_DEN_genTracks"].evaluate(abseta_flat, pt_flat, "systdown")

    reco_nom_flat = ceval_HighPt["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evaluate(abseta_flat, pt_flat, "nominal")
    reco_up_flat = ceval_HighPt["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evaluate(abseta_flat, pt_flat, "systup")
    reco_down_flat = ceval_HighPt["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evaluate(abseta_flat, pt_flat, "systdown")

    ID_nom_flat = ceval_HighPt["NUM_HighPtID_DEN_GlobalMuonProbes"].evaluate(abseta_flat, pt_flat, "nominal")
    ID_up_flat = ceval_HighPt["NUM_HighPtID_DEN_GlobalMuonProbes"].evaluate(abseta_flat, pt_flat, "systup")
    ID_down_flat = ceval_HighPt["NUM_HighPtID_DEN_GlobalMuonProbes"].evaluate(abseta_flat, pt_flat, "systdown")
    
    ISO_nom_flat = ceval_HighPt["NUM_probe_LooseRelTkIso_DEN_HighPtProbes"].evaluate(abseta_flat, pt_flat, "nominal") 
    ISO_up_flat = ceval_HighPt["NUM_probe_LooseRelTkIso_DEN_HighPtProbes"].evaluate(abseta_flat, pt_flat, "systup") 
    ISO_down_flat = ceval_HighPt["NUM_probe_LooseRelTkIso_DEN_HighPtProbes"].evaluate(abseta_flat, pt_flat, "systdown") 

    # unflatten arrays

    tracking_nom = ak.unflatten(tracking_nom_flat, ak.num(pt))
    tracking_up = ak.unflatten(tracking_up_flat, ak.num(pt))
    tracking_down = ak.unflatten(tracking_down_flat, ak.num(pt))

    reco_nom = ak.unflatten(reco_nom_flat, ak.num(pt))
    reco_up = ak.unflatten(reco_up_flat, ak.num(pt))
    reco_down = ak.unflatten(reco_down_flat, ak.num(pt))

    ID_nom = ak.unflatten(ID_nom_flat, ak.num(pt))
    ID_up = ak.unflatten(ID_up_flat, ak.num(pt))
    ID_down = ak.unflatten(ID_down_flat, ak.num(pt))

    ISO_nom = ak.unflatten(ISO_nom_flat, ak.num(pt))
    ISO_up = ak.unflatten(ISO_up_flat, ak.num(pt))
    ISO_down = ak.unflatten(ISO_down_flat, ak.num(pt))

    # attach SFs to muons 
    muons['SF_muon_nom'] = tracking_nom * reco_nom * ID_nom * ISO_nom
    muons['SF_muon_up'] = tracking_up * reco_up * ID_up * ISO_up
    muons['SF_muon_down'] = tracking_down * reco_down * ID_down * ISO_down

    muons['SF_ele_nom'] = ak.ones_like(pt)
    muons['SF_ele_up'] = ak.ones_like(pt)
    muons['SF_ele_down'] = ak.ones_like(pt)    


def ApplyMuonPtCorr(muons, year, isData):
    '''
    For muons with pt<120, Rochester Corrections are applied 
    For muons with pt>120, the tunep correction is applied
    '''

    pt = muons.pt
    pt_flat = ak.flatten(pt)
    pt_mask = ak.flatten(pt >= 120)
    tunep_factor = ak.flatten(muons.tunepRelPt)

    rocco_pt = ak.flatten(ApplyRochesterCorrections(muons, year, isData))
    tunep_pt = pt_flat * tunep_factor

    pt_corr_flat = ak.where(pt_mask, rocco_pt, tunep_pt) #if pt mask, pt_corr=tunep_pt, else pt_corr=rocco_pt
    pt_corr = ak.unflatten(pt_corr_flat, ak.num(pt))

    return pt_corr


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


###### Muon Rochester corrections
################################################################
# https://gitlab.cern.ch/akhukhun/roccor
# https://github.com/CoffeaTeam/coffea/blob/master/coffea/lookup_tools/rochester_lookup.py
def ApplyRochesterCorrections(mu, year, isData):
    rocco_tag = None

    rocco_year_map = {
        '2016APV': '2016aUL',
        '2016': '2016bUL',
        '2017': '2017UL',
        '2018': '2018UL'
    }

    rocco_tag = rocco_year_map[year]
    rochester_data = txt_converters.convert_rochester_file(topcoffea_path(f"data/MuonScale/RoccoR{rocco_tag}.txt"), loaduncs=True)
    rochester = rochester_lookup.rochester_lookup(rochester_data)
    if not isData:
        hasgen = ~np.isnan(ak.fill_none(mu.matched_gen.pt, np.nan))
        mc_rand = np.random.rand(*ak.to_numpy(ak.flatten(mu.pt)).shape)
        mc_rand = ak.unflatten(mc_rand, ak.num(mu.pt, axis=1))
        corrections = np.array(ak.flatten(ak.ones_like(mu.pt)))
        mc_kspread = rochester.kSpreadMC(
            mu.charge[hasgen],mu.pt[hasgen],
            mu.eta[hasgen],
            mu.phi[hasgen],
            mu.matched_gen.pt[hasgen]
        )
        mc_ksmear = rochester.kSmearMC(
            mu.charge[~hasgen],
            mu.pt[~hasgen],
            mu.eta[~hasgen],
            mu.phi[~hasgen],
            mu.nTrackerLayers[~hasgen],
            mc_rand[~hasgen]
        )
        hasgen_flat = np.array(ak.flatten(hasgen))
        corrections[hasgen_flat] = np.array(ak.flatten(mc_kspread))
        corrections[~hasgen_flat] = np.array(ak.flatten(mc_ksmear))
        corrections = ak.unflatten(corrections, ak.num(mu.pt, axis=1))
    else:
        corrections = rochester.kScaleDT(mu.charge, mu.pt, mu.eta, mu.phi)

    return (mu.pt * corrections)
