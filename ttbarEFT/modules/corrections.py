import numpy as np
import pickle
import gzip
import awkward as ak
import re
import yaml
import uproot

import correctionlib

from coffea import lookup_tools
from coffea.lookup_tools import txt_converters, rochester_lookup

from topcoffea.modules.paths import topcoffea_path
from ttbarEFT.modules.paths import ttbarEFT_path
from ttbarEFT.modules.CorrectedJetsFactory import CorrectedJetsFactory #, get_jec_uncertainty_label
from ttbarEFT.modules.CorrectedMETFactory import CorrectedMETFactory
from ttbarEFT.modules.JECStack import JECStack


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

goldenJSON_map = {
    "2016APV": "Collisions16_UltraLegacy_goldenJSON",
    "2016": "Collisions16_UltraLegacy_goldenJSON",
    "2017": "Collisions17_UltraLegacy_goldenJSON",
    "2018": "Collisions18_UltraLegacy_goldenJSON",
    "2022": "Collisions2022_355100_357900_eraBCD_GoldenJson",
    "2022EE": "Collisions2022_359022_362760_eraEFG_GoldenJson",
    "2023": "Collisions2023_366403_369802_eraBC_GoldenJson",
    "2023BPix": "Collisions2023_369803_370790_eraD_GoldenJson",
}


def GetPUSF(nTrueInt, year, var='nominal'):
    year = str(year)
    if year not in clib_year_map.keys():
        raise Exception(f"Error: Unknown year \"{year}\".")

    clib_year = clib_year_map[year]
    json_path = ttbarEFT_path(f"data/POG/LUM/{clib_year}/puWeights.json.gz")
    ceval = correctionlib.CorrectionSet.from_file(json_path)

    pucorr_tag = goldenJSON_map[year]
    pu_corr = ceval[pucorr_tag].evaluate(nTrueInt, var)
    return pu_corr


def AttachScaleWeights(events):
    """
    Dynamically retrieves scale weights from LHEScaleWeight based on its __doc__.

    LHE scale variation weights (w_var / w_nominal)
    """

    # Check if LHEScaleWeight exists in the event
    if events.LHEScaleWeight is None:
        raise Exception('LHEScaleWeight not found!')

    # Get the LHEScaleWeight documentation
    scale_weight_doc = events.LHEScaleWeight.__doc__

    does_doc_exist = True
    if not scale_weight_doc:
        #raise Exception('LHEScaleWeight.__doc__ is empty or not available!')
        does_doc_exist = False

    # Define the mapping we are looking for the three scenarios
    scenarios_map = {
        # Scenario 1: renscfact and facscfact for 9 weights
        "renscfact": {
            "scale_map": {
                'renscfact=0.5d0 facscfact=0.5d0': 'renormfactDown',
                'renscfact=0.5d0 facscfact=1d0': 'renormDown',
                'renscfact=0.5d0 facscfact=2d0': 'renormDown_factUp',
                'renscfact=1d0 facscfact=0.5d0': 'factDown',
                #'renscfact=1d0 facscfact=1d0': 'nominal',  # Explicitly handle nominal case in processor
                'renscfact=1d0 facscfact=2d0': 'factUp',
                'renscfact=2d0 facscfact=0.5d0': 'renormUp_factDown',
                'renscfact=2d0 facscfact=1d0': 'renormUp',
                'renscfact=2d0 facscfact=2d0': 'renormfactUp'
            },
            "re_pattern": r'\[(\d+)\] is renscfact=(\d+\.?\d*)d0 facscfact=(\d+\.?\d*)d0',
            "key": lambda match: f'renscfact={match[1]}d0 facscfact={match[2]}d0'
        },
        # Scenario 2: MUF and MUR for 9 weights
        "MUF9": {
            "scale_map": {
                'MUF="0.5" MUR="0.5"': 'renormfactDown',
                'MUF="1.0" MUR="0.5"': 'renormDown',
                'MUF="2.0" MUR="0.5"': 'renormDown_factUp',
                'MUF="0.5" MUR="1.0"': 'factDown',
                #'MUF="1.0" MUR="1.0"': 'nominal',  #  Explicitly handle nominal case in processor
                'MUF="2.0" MUR="1.0"': 'factUp',
                'MUF="0.5" MUR="2.0"': 'renormUp_factDown',
                'MUF="1.0" MUR="2.0"': 'renormUp',
                'MUF="2.0" MUR="2.0"': 'renormfactUp'
            },
            "re_pattern": r'\[(\d+)\] is MUF="(\d+\.?\d*)" MUR="(\d+\.?\d*)"',
            "key": lambda match: f'MUF="{match[1]}" MUR="{match[2]}"'
        },
        # Scenario 3: MUF and MUR for 8 weights
        "MUF8": {
            "scale_map": {
                'MUF="0.5" MUR="0.5"': 'renormfactDown',
                'MUF="1.0" MUR="0.5"': 'renormDown',
                'MUF="2.0" MUR="0.5"': 'renormDown_factUp',
                'MUF="0.5" MUR="1.0"': 'factDown',
                'MUF="2.0" MUR="1.0"': 'factUp',
                'MUF="0.5" MUR="2.0"': 'renormUp_factDown',
                'MUF="1.0" MUR="2.0"': 'renormUp',
                'MUF="2.0" MUR="2.0"': 'renormfactUp'
            },
            "re_pattern": r'\[(\d+)\] is MUF="(\d+\.?\d*)" MUR="(\d+\.?\d*)"',
            "key": lambda match: f'MUF="{match[1]}" MUR="{match[2]}"'
        }
    }

    # Determine the number of weights available
    len_of_wgts = ak.count(events.LHEScaleWeight, axis=-1)
    all_len_9_or_0_bool = ak.all((len_of_wgts == 9) | (len_of_wgts == 0))
    all_len_8_or_0_bool = ak.all((len_of_wgts == 8) | (len_of_wgts == 0))
    scale_weights = None
    scale_map = None
    matches = None
    scenario = None

    # Choose between the different cases based on the number of weights and the doc string
    if all_len_9_or_0_bool:
        if "renscfact" in scale_weight_doc:
            scenario = "renscfact"  # Scenario 1: renscfact/facscfact
        elif "MUF" in scale_weight_doc:
            scenario = "MUF9"       # Scenario 2: MUF/MUR with 9 weights
    elif all_len_8_or_0_bool:
        scenario = "MUF8"           # Scenario 3: MUF/MUR with 8 weights
    else:
        raise Exception("Unknown weight type")

    scale_weights = ak.fill_none(ak.pad_none(events.LHEScaleWeight, 9 if (scenario == "MUF9" or scenario == "renscfact" or scenario is None) else 8), 1)
    # Dictionary to hold the index of each variation
    scale_indices = {}

    if scenario is not None:
        matches = re.findall(scenarios_map[scenario]["re_pattern"], scale_weight_doc)
        scale_map = scenarios_map[scenario]["scale_map"]
        key = scenarios_map[scenario]["key"]

        # Parse the matches and build the scale indices dictionary
        for match in matches:
            index = int(match[0])  # Extract the index from the regex match
            key_str = key(match)  # Dynamically get the key string based on the case

            if key_str in scale_map:
                scale_indices[scale_map[key_str]] = index

        required_keys = list(scale_map.values())

        # Check if all needed weights were found
        if not all(key in scale_indices for key in required_keys):
            missing_keys = [key for key in required_keys if key not in scale_indices]
            raise Exception('Not all scale weight variations found in LHEScaleWeight.__doc__!')

    else:
        #This part of the code assumes that every entry is unit when LHEScaleWeight is not actually filled
        dummy_keys = list(scenarios_map["MUF9"]["scale_map"].values())
        for id_key, dummy_key in enumerate(dummy_keys):
            scale_indices[dummy_key] = id_key

    # Assign the weights from the event to the respective fields dynamically using a loop
    for key in scale_indices:
        events[key] = scale_weights[:, scale_indices[key]]



############################
######### btag SF  #########
############################

def GetBtagEff(year, jets, wp='medium'):
    # similar to GetMCeffFunc in topeft.modules.corrections
    if year not in clib_year_map.keys():
        raise Exception(f"Error: Unknown year \"{year}\".")

    pathToBtagMCeff = ttbarEFT_path(f"data/btagSF/UL/btagMCeff_{year}.pkl.gz")
    hists = {}
    with gzip.open(pathToBtagMCeff) as fin:
        hists = pickle.load(fin)['btag']

    h = hists['jetptetaflav']
    hnum = h[{'WP': wp}]
    hden = h[{'WP': 'all'}]

    eff = hnum/hden
    eff_lookup = lookup_tools.dense_lookup.dense_lookup(
        eff.values(), 
        [ax.edges for ax in eff.axes]
    )

    # this order must match the order of ax in eff, which is based on btagMCeff_processor order when the hist is initialized
    return eff_lookup(jets.hadronFlavour, jets.pt, np.abs(jets.eta))


def GetBtagEffLookup(year, wp='medium'):

    if year not in clib_year_map.keys(): 
        raise Exception(f"Error: Unknown year \"{year}\".")

    pathToBtagMCeff = ttbarEFT_path(f'data/btagSF/UL/btagMCeff_{year}.pkl.gz')
    
    hists = {}
    with gzip.open(pathToBtagMCeff) as fin:
        hists = pickle.load(fin)['btag']

    h = hists['jetptetaflav']
    hnum = h[{'WP': wp}]
    hden = h[{'WP': 'all'}]

    eff = hnum/hden
    eff_lookup = lookup_tools.dense_lookup.dense_lookup(
        eff.values(), 
        [ax.edges for ax in eff.axes]
    )

    def evaluate_eff(jets):
        return eff_lookup(jets.hadronFlavour, jets.pt, np.abs(jets.eta))

    return evaluate_eff


def GetBtagSF(jet_collection,wp,year,method,syst):
    """
    Get btag SF from central correctionlib json
    - similar to topcoffea.modules.corrections.btag_sf_eval()
    - usage: btag_sfL_up   = tc_cor.btag_sf_eval(jets_flav, "L",sys_year,f"deepJet_{dJ_tag}",f"up_{corrtype}")
    - usage: weights_obj_base_for_kinematic_syst.add(f"btagSF{b_syst}", events.nom, btag_w_up, btag_w_down)
    """

    clib_year = clib_year_map[year]
    fname = ttbarEFT_path(f"data/POG/BTV/{clib_year}/btagging.json.gz")
    ceval = correctionlib.CorrectionSet.from_file(fname)

    # Flatten the input (until correctionlib handles jagged data natively)
    abseta_flat = ak.flatten(abs(jet_collection.eta))
    pt_flat = ak.flatten(jet_collection.pt)
    flav_flat = ak.flatten(jet_collection.hadronFlavour)

    # For now, cap all pt at 1000 https://cms-talk.web.cern.ch/t/question-about-evaluating-sfs-with-correctionlib/31763
    pt_flat = ak.where(pt_flat>1000.0,1000.0,pt_flat)

    # Evaluate the SF
    sf_flat = ceval[method].evaluate(syst,wp,flav_flat,abseta_flat,pt_flat)
    sf = ak.unflatten(sf_flat,ak.num(jet_collection.pt))

    return sf


def GetBtagSFLookup(wp,year,method):
    """
    Get btag SF from central correctionlib json
    - similar to topcoffea.modules.corrections.btag_sf_eval()
    - usage: btag_sfL_up   = tc_cor.btag_sf_eval(jets_flav, "L",sys_year,f"deepJet_{dJ_tag}",f"up_{corrtype}")
    - usage: weights_obj_base_for_kinematic_syst.add(f"btagSF{b_syst}", events.nom, btag_w_up, btag_w_down)
    """

    clib_year = clib_year_map[year]
    fname = ttbarEFT_path(f"data/POG/BTV/{clib_year}/btagging.json.gz")
    ceval = correctionlib.CorrectionSet.from_file(fname)[method]

    def evaluate_sf(jet_collection, syst):
        # Flatten the input (until correctionlib handles jagged data natively)
        abseta_flat = ak.flatten(abs(jet_collection.eta))
        pt_flat = ak.flatten(jet_collection.pt)
        flav_flat = ak.flatten(jet_collection.hadronFlavour)

        # For now, cap all pt at 1000 https://cms-talk.web.cern.ch/t/question-about-evaluating-sfs-with-correctionlib/31763
        pt_flat = ak.where(pt_flat>1000.0,1000.0,pt_flat)

        # Evaluate the SF
        sf_flat = ceval.evaluate(syst,wp,flav_flat,abseta_flat,pt_flat)
        sf = ak.unflatten(sf_flat,ak.num(jet_collection.pt))

        return sf

    return evaluate_sf


def GetBtag_method1a_wgt_singlewp(eff,sf,passes_tag):
    """
    Evaluate btag method 1a weight for a single WP (https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods)
    - Takes as input a given array of eff and sf and a mask for whether or not the events pass a tag
    - Returns P(DATA)/P(MC)
    - Where P(MC) = Product over tagged (eff) * Product over not tagged (1-eff)
    - Where P(DATA) = Product over tagged (eff*sf) * Product over not tagged (1-eff*sf)
    """

    p_mc = ak.prod(eff[passes_tag],axis=-1) * ak.prod(1-eff[~passes_tag],axis=-1)
    p_data = ak.prod(eff[passes_tag]*sf[passes_tag],axis=-1) * ak.prod(1-eff[~passes_tag]*sf[~passes_tag],axis=-1)
    wgt = p_data/p_mc

    return wgt



###################################
######### Jet Corrections #########
###################################

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


def GetJetPuIDSF(year, jets, var='nom', wp='L'):

    year = str(year)
    if year not in clib_year_map.keys():
        raise Exception(f"Error: Unknown year \"{year}\".")

    clib_year = clib_year_map[year]
    json_path = ttbarEFT_path(f"data/POG/JME/{clib_year}/jmar.json.gz")
    ceval = correctionlib.CorrectionSet.from_file(json_path)
    key = "PUJetID_eff"

    eta_flat = ak.flatten(jets.eta)
    pt_flat = ak.flatten(jets.pt)
    mask_flat = pt_flat < 50 

    pt_flat_for_eval = ak.where(mask_flat, pt_flat, 49.9)

    weights_raw = ceval[key].evaluate(eta_flat, pt_flat_for_eval, var, wp)
    jetPuID_flat = ak.where(mask_flat, weights_raw, 1.0)
    jetPuID_unflat = ak.unflatten(jetPuID_flat, ak.num(jets.pt))

    return ak.prod(jetPuID_unflat, axis=1)



def get_jerc_keys(year, isdata, era=None):

    #JERC dictionary for various keys
    with open(ttbarEFT_path("modules/jerc_dict.yaml"), "r") as f:
        jerc_dict = yaml.safe_load(f)

    # Jet Algorithm
    if year.startswith("202"):
        jet_algo = 'AK4PFPuppi'
        jec_levels = jerc_dict[year]['jec_levels']
    else:
        jet_algo = 'AK4PFchs'
        jec_levels = []            # no JEC corrections for Run2, already applied in nanoAODv9
    

    # jerc keys and junc types
    if not isdata:
        jec_key    = jerc_dict[year]['jec_mc']
        jer_key    = jerc_dict[year]['jer']
        junc_types = jerc_dict[year]['junc']
    else:
        # if year in ['2016']: #,'2022','2023BPix'
        #     jec_key = jerc_dict[year]['jec_data']
        # else:
        #     jec_key = jerc_dict[year]['jec_data'][era]
        jec_key     = None
        jer_key     = None
        junc_types  = None

    return jet_algo, jec_key, jec_levels, jer_key, junc_types


def get_supported_jet_systematics(year, isData=False, era=None):
    if isData:
        return []

    systs = [f"JER_{year}Up", f"JER_{year}Down"]
    for base in get_supported_jes_bases(year, isData=isData, era=era):
        systs.append(f"JES_{base}Up")
        systs.append(f"JES_{base}Down")
    return systs



def get_supported_jes_bases(year, isData=False, era=None):
    if isData:
        return []

    jet_algo, jec_tag, _, _, junc_types = get_jerc_keys(year, isData, era)

    bases = []
    for junc_type in junc_types:
        full_unc_name = f"{jec_tag}_{junc_type}_{jet_algo}"
        bases.append(get_jec_uncertainty_label(full_unc_name, jec_tag, jet_algo))

    if len(set(bases)) != len(bases):
        duplicates = sorted({x for x in bases if bases.count(x) > 1})
        raise RuntimeError(
            f"Collision in JES base labels for year {year}: {duplicates}"
        )

    return bases


def get_jec_uncertainty_label(junc_name, jec_tag, jet_algo):
    """
    Extract the uncertainty label from a full correction name:
      {jec_tag}_{junc_type}_{jet_algo}
    """
    prefix = f"{jec_tag}_"
    suffix = f"_{jet_algo}"
    if not junc_name.startswith(prefix):
        raise ValueError(
            f'Uncertainty name "{junc_name}" does not start with expected prefix "{prefix}".'
        )
    if not junc_name.endswith(suffix):
        raise ValueError(
            f'Uncertainty name "{junc_name}" does not end with expected suffix "{suffix}".'
        )
    junc_type = junc_name[len(prefix) : -len(suffix)]
    if not junc_type:
        raise ValueError(f'Failed to extract junc_type from uncertainty name "{junc_name}".')
    return _canonicalize_junc_type_label(junc_type, jec_tag)


def _is_run2_jec_tag(jec_tag):
    return "UL" in jec_tag


def _canonicalize_junc_type_label(junc_type, jec_tag):
    # Preserve existing Run2 naming: regrouped UL sources drop the "Regrouped_" prefix.
    if _is_run2_jec_tag(jec_tag) and junc_type.startswith("Regrouped_"):
        return junc_type.replace("Regrouped_", "", 1)
    return junc_type


####### JEC
##############################################
# JER: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
# JES: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC

def ApplyJetCorrections(year, corr_type, isData, era, useclib=True, savelevels=False):

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
    
    usejecstack = not useclib

    if year not in clib_year_map.keys():
        raise Exception(f"Error: Unknown year \"{year}\".")

    jec_year = clib_year_map[year]

    jet_algo, jec_tag, jec_levels, jer_tag, junc_types = get_jerc_keys(year, isData, era)
    json_path = ttbarEFT_path(f"data/POG/JME/{jec_year}/jet_jerc.json.gz")

    # Create JECStack for clib scenario
    jec_stack = JECStack(
        jec_tag=jec_tag,
        jec_levels=[],      # should be jec_levels for Run3 or if JEC corrections need to be applied
        jer_tag=jer_tag,
        jet_algo=jet_algo,
        junc_types=junc_types,
        json_path=json_path,
        use_clib=useclib,
        savecorr=savelevels
    )

    # Name map for jet or MET corrections
    name_map = {
        'JetPt':    'pt',
        'JetMass':  'mass',
        'JetEta':   'eta',
        'JetPhi':   'phi',
        'JetA':     'area',
        'ptGenJet': 'pt_gen',
        'ptRaw':    'pt',       # should be pt_raw for Run3 or if JEC corrections need to be added 
        'massRaw':  'mass',     # should be mass_raw for Run3 or if JEC corrections need to be added 
        'Rho':      'Rho',
        'METpt':    'pt',
        'METphi':   'phi',
        'UnClusteredEnergyDeltaX': 'MetUnclustEnUpDeltaX',
        'UnClusteredEnergyDeltaY': 'MetUnclustEnUpDeltaY'
    }

    MET_name_map = {
        'JetPt':    'pt',
        'JetPhi':   'phi',
        'ptRaw':    'pt_orig',       
        'METpt':    'pt',
        'METphi':   'phi',
        'UnClusteredEnergyDeltaX': 'MetUnclustEnUpDeltaX',
        'UnClusteredEnergyDeltaY': 'MetUnclustEnUpDeltaY'
    }

    # Return appropriate factory based on correction type
    if corr_type == 'met':
        return CorrectedMETFactory.CorrectedMETFactory(name_map)
    elif corr_type == 'jets':
        return CorrectedJetsFactory(name_map, jec_stack, allowed_variations=None)
    else: 
        raise ValueError(f"Unknown correction type \"{corr_type}\".")


def ApplyJetSystematics(year,cleanedJets,syst_var):

    # getattr(cleanedJets, f"JER.up.pt")
    if (syst_var == 'nominal'):
        return cleanedJets

    if 'Up' in syst_var: 
        suffix = 'up'
        if syst_var.startswith(f'JER_{year}'):
            base_field = "JER"
        elif syst_var == "JESUp":
            base_field = "JES_jes"
        else: 
            base_field = syst_var[:-2]
            
    elif 'Down' in syst_var:
            suffix = 'down'
            if syst_var.startswith(f'JER_{year}'):
                base_field = "JER"
            elif syst_var == "JESDown":
                base_field = "JES_jes"
            else:
                base_field = syst_var[:-4]
    else: 
        raise ValueError(f"Unknown variation {syst_var}.")
        
    if base_field not in cleanedJets.fields:
        raise ValueError(f"Field {base_field} not found in jets for variation {syst_var}")
    
    var_record = getattr(cleanedJets[base_field], suffix)
    new_jets = ak.with_field(cleanedJets, var_record["pt"], "pt")
    new_jets = ak.with_field(new_jets, var_record["mass"], "mass")
    
    return new_jets

    # # Save `2016APV` as `2016APV` but look up `2016` corrections (no separate APV corrections available)
    # elif ('Up' in syst_var and syst_var[:-2].replace('APV', '') in cleanedJets.fields):
    #     return cleanedJets[syst_var.replace('Up', '').replace("Pile", "PileUp").replace('APV', '')].up
    # elif ('Down' in syst_var and syst_var[:-4].replace('APV', '') in cleanedJets.fields):
    #     return cleanedJets[syst_var.replace('Down', '').replace('APV', '')].down
    # else:
    #     raise Exception(f"Error: Unknown variation \"{syst_var}\".")


########################################
######### Lepton Scale Factors #########
########################################

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
    electrons['SF_ele_up'] = reco_up * HEEPSF_nom * HEEPSF_up 
    electrons['SF_ele_down'] = reco_down * HEEPSF_nom * HEEPSF_down 

    # electrons['SF_eleRECO_nom'] = reco_nom 
    # electrons['SF_eleRECO_up'] = reco_up 
    # electrons['SF_eleRECO_down'] = reco_down 

    # electrons['SF_eleHEEP_nom'] = HEEPSF_nom
    # electrons['SF_eleHEEP_up'] = HEEPSF_nom * HEEPSF_up 
    # electrons['SF_eleHEEP_down'] = HEEPSF_nom * HEEPSF_down 

    electrons['SF_muonID_nom'] = ak.ones_like(reco_nom)
    electrons['SF_muonID_up'] = ak.ones_like(reco_nom)
    electrons['SF_muonID_down'] = ak.ones_like(reco_nom)

    electrons['SF_muonISO_nom'] = ak.ones_like(reco_nom)
    electrons['SF_muonISO_up'] = ak.ones_like(reco_nom)
    electrons['SF_muonISO_down'] = ak.ones_like(reco_nom)


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

    muonRECO_year_map = {
        "2016APV": "2016_preVFP",
        "2016preVFP": "2016_preVFP",
        "2016": "2016",
        "2017": "2017",
        "2018": "2018",
    }

    # initialize muon variables 
    abs_eta = np.abs(muons.eta) #For run3 abs(eta) should be changed to signed eta
    pt = muons.pt 
    phi = muons.phi 

    abseta_flat = ak.flatten(abs_eta)
    pt_flat = ak.flatten(pt)
    p_flat = np.multiply(pt_flat, np.cosh(abseta_flat))
    phi_flat = ak.flatten(phi)

    clib_year = clib_year_map[year]
    muonRECO_year = muonRECO_year_map[year]
    # json_path_HighPt = ttbarEFT_path(f"data/POG/MUO/{clib_year}/muon_HighPt.json.gz")
    json_path_z = ttbarEFT_path(f"data/POG/MUO/{clib_year}/muon_Z.json.gz")
    ceval_z = correctionlib.CorrectionSet.from_file(json_path_z)
    ceval_RECO = correctionlib.CorrectionSet.from_file(ttbarEFT_path(f"data/POG/MUO/{clib_year}/ScaleFactors_Muon_highPt_RECO_{muonRECO_year}_schemaV2.json"))
    ceval_IDISO = correctionlib.CorrectionSet.from_file(ttbarEFT_path(f"data/POG/MUO/{clib_year}/ScaleFactors_Muon_highPt_IDISO_{clib_year}_schemaV2.json"))
    
    # pt_mask = ak.flatten(pt >= 50) # the lowest pT bin in muon_HighPt is 50 #
    # evaluate SFs
    tracking_nom_flat = ceval_z["NUM_TrackerMuons_DEN_genTracks"].evaluate(abseta_flat, pt_flat, "nominal")
    tracking_up_flat = ceval_z["NUM_TrackerMuons_DEN_genTracks"].evaluate(abseta_flat, pt_flat, "systup")
    tracking_down_flat = ceval_z["NUM_TrackerMuons_DEN_genTracks"].evaluate(abseta_flat, pt_flat, "systdown")

    # high pT RECO SF are binned in eta and momentum, not eta and pT
    reco_nom_flat = ceval_RECO["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evaluate(abseta_flat, p_flat, "nominal")
    reco_up_flat = ceval_RECO["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evaluate(abseta_flat, p_flat, "systup")
    reco_down_flat = ceval_RECO["NUM_GlobalMuons_DEN_TrackerMuonProbes"].evaluate(abseta_flat, p_flat, "systdown")

    ID_nom_flat = ceval_IDISO["NUM_HighPtID_DEN_GlobalMuonProbes"].evaluate(abseta_flat, pt_flat, "nominal")
    ID_up_flat = ceval_IDISO["NUM_HighPtID_DEN_GlobalMuonProbes"].evaluate(abseta_flat, pt_flat, "systup")
    ID_down_flat = ceval_IDISO["NUM_HighPtID_DEN_GlobalMuonProbes"].evaluate(abseta_flat, pt_flat, "systdown")
    
    ISO_nom_flat = ceval_IDISO["NUM_probe_LooseRelTkIso_DEN_HighPtProbes"].evaluate(abseta_flat, pt_flat, "nominal") 
    ISO_up_flat = ceval_IDISO["NUM_probe_LooseRelTkIso_DEN_HighPtProbes"].evaluate(abseta_flat, pt_flat, "systup") 
    ISO_down_flat = ceval_IDISO["NUM_probe_LooseRelTkIso_DEN_HighPtProbes"].evaluate(abseta_flat, pt_flat, "systdown") 

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
    muons['SF_muonID_nom'] = tracking_nom * ID_nom * reco_nom
    muons['SF_muonID_up'] = tracking_up * ID_up * reco_up
    muons['SF_muonID_down'] = tracking_down * ID_down * reco_down

    muons['SF_muonISO_nom'] = ISO_nom
    muons['SF_muonISO_up'] = ISO_up
    muons['SF_muonISO_down'] = ISO_down

    muons['SF_ele_nom'] = ak.ones_like(pt)
    muons['SF_ele_up'] = ak.ones_like(pt)
    muons['SF_ele_down'] = ak.ones_like(pt)   

    # muons['SF_eleRECO_nom'] = ak.ones_like(pt)
    # muons['SF_eleRECO_up'] = ak.ones_like(pt)
    # muons['SF_eleRECO_down'] = ak.ones_like(pt)

    # muons['SF_eleHEEP_nom'] = ak.ones_like(pt)
    # muons['SF_eleHEEP_up'] = ak.ones_like(pt)
    # muons['SF_eleHEEP_down'] = ak.ones_like(pt)


def Get_ElecIDSF(events):

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]
 
    calc_nom = l0.SF_ele_nom * l1.SF_ele_nom
    calc_up = l0.SF_ele_up * l1.SF_ele_up
    calc_down = l0.SF_ele_down * l1.SF_ele_down

    return ak.fill_none(calc_nom, 1.0), ak.fill_none(calc_up, 1.0), ak.fill_none(calc_down, 1.0)


def Get_ElecRECOSF(events):

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]

    calc_nom = l0.SF_eleRECO_nom * l1.SF_eleRECO_nom
    calc_up = l0.SF_eleRECO_up * l1.SF_eleRECO_up
    calc_down = l0.SF_eleRECO_down * l1.SF_eleRECO_down

    return ak.fill_none(calc_nom, 1.0), ak.fill_none(calc_up, 1.0), ak.fill_none(calc_down, 1.0)


def Get_MuonIDSF(events):

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]

    calc_nom = l0.SF_muonID_nom * l1.SF_muonID_nom
    calc_up = l0.SF_muonID_up * l1.SF_muonID_up
    calc_down = l0.SF_muonID_down * l1.SF_muonID_down

    return ak.fill_none(calc_nom, 1.0), ak.fill_none(calc_up, 1.0), ak.fill_none(calc_down, 1.0)


def Get_MuonISOSF(events):

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]

    calc_nom = l0.SF_muonISO_nom * l1.SF_muonISO_nom
    calc_up = l0.SF_muonISO_up * l1.SF_muonISO_up
    calc_down = l0.SF_muonISO_down * l1.SF_muonISO_down

    return ak.fill_none(calc_nom, 1.0), ak.fill_none(calc_up, 1.0), ak.fill_none(calc_down, 1.0)


def GetLepSF(events, lep_cat):

    # Select events based on lepton category
    mask = events[f'is_{lep_cat}']

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]

    def get_flavor_sf(lep, var):
        return ak.where(abs(lep.pdgId) == 11, 
                        getattr(lep, f"SF_ele_{var}"), 
                        getattr(lep, f"SF_muon_{var}"))

    if lep_cat == 'ee': 
        calc_nom = l0.SF_ele_nom * l1.SF_ele_nom
        calc_up = l0.SF_ele_up * l1.SF_ele_up
        calc_down = l0.SF_ele_down * l1.SF_ele_down

    elif lep_cat == 'mm': 
        calc_nom = l0.SF_muon_nom * l1.SF_muon_nom
        calc_up = l0.SF_muon_up * l1.SF_muon_up
        calc_down = l0.SF_muon_down * l1.SF_muon_down

    elif lep_cat == 'em':
        calc_nom = get_flavor_sf(l0, 'nom') * get_flavor_sf(l1, 'nom')
        calc_up = get_flavor_sf(l0, 'up') * get_flavor_sf(l1, 'up')
        calc_down = get_flavor_sf(l0, 'down') * get_flavor_sf(l1, 'down')

    nom = ak.where(mask, calc_nom, 1.0)
    up = ak.where(mask, calc_up, 1.0)
    down = ak.where(mask, calc_down, 1.0)

    return nom, up, down


####################################
######### Muon Corrections #########
####################################

def ApplyMuonPtCorr(muons, year, isData):
    '''
    For muons with pt<120, Rochester Corrections are applied 
    For all muons, the tunep correction is applied
    Returns corrected pt
    '''

    pt = muons.pt
    pt_flat = ak.flatten(pt)
    pt_mask = ak.flatten(pt >= 120)
    tunep_factor = ak.flatten(muons.tunepRelPt)

    rocco_pt = ak.flatten(ApplyRochesterCorrections(muons, year, isData)) #output of rocco is correction*orig_pt
    tunep_pt = tunep_factor * pt_flat   # calculate correction*orig_pt

    pt_corr_flat = ak.where(pt_mask, tunep_pt, rocco_pt*tunep_factor) #if pt>120, pt_=tunep_pt, else pt=rocco_pt*tunep_factor
    pt_corr = ak.unflatten(pt_corr_flat, ak.num(pt))

    return pt_corr


def ApplyRochesterCorrections(mu, year, isData):
    # Muon Rochester corrections
    # https://gitlab.cern.ch/akhukhun/roccor
    # https://github.com/CoffeaTeam/coffea/blob/master/coffea/lookup_tools/rochester_lookup.py

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


#####################################
######## Trigger Efficienies ########
#####################################

def AttachElecTrigEff_nouncert(electrons, year):

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


    eta_flat = ak.flatten(eta)
    Et_flat = ak.flatten(electrons.pt * (electrons.scEtOverPt+1))

    ElecSF_flat = ak.where(
                eta_barrel_mask, 
                SFevaluator[f"ElecSF_{year}_barrel"](Et_flat),
                SFevaluator[f"ElecSF_{year}_endcap"](Et_flat)
            )

    electrons['trig_eff_ele_nom'] = ak.unflatten(ElecSF_flat, ak.num(eta))

    # fill muon trig efficiencies as 1.0
    electrons['trig_MCeff_mu_nom'] = ak.ones_like(electrons.pt)
    electrons['trig_DATAeff_mu_nom'] = ak.ones_like(electrons.pt)
    electrons['trig_MCeff_mu_up'] = ak.ones_like(electrons.pt)
    electrons['trig_DATAeff_mu_up'] = ak.ones_like(electrons.pt)
    electrons['trig_MCeff_mu_down'] = ak.ones_like(electrons.pt)
    electrons['trig_DATAeff_mu_down'] = ak.ones_like(electrons.pt)


# new function that includes uncertainties
def AttachElecTrigEff(electrons, year):
    eff_path = ttbarEFT_path('data/leptonSF/elec/DiEleCaloIdLMWPMS2_HEEPeff.root')
    
    # Map input year to ROOT internal tags
    year_map = {
        "2016APV": "UL2016preVFP",
        "2016":    "UL2016postVFP",
        "2017":    "UL2017",
        "2018":    "UL2018",
    }
    tag = year_map[year]

    # Extract values and variances manually
    with uproot.open(eff_path) as f:
        # Barrel
        h_b = f[f"{tag}_Barrel_Et"]
        edges_b = h_b.axes[0].edges()
        val_b = h_b.values()
        err_b = np.sqrt(h_b.variances()) #Convert variance to sigma
        
        # Endcap
        h_e = f[f"{tag}_Endcaps_Et"]
        edges_e = h_e.axes[0].edges()
        val_e = h_e.values()
        err_e = np.sqrt(h_e.variances())

    # Create the dense_lookup objects for nominal and errors
    lookup_sf_b = lookup_tools.dense_lookup.dense_lookup(val_b, edges_b)
    lookup_sf_e = lookup_tools.dense_lookup.dense_lookup(val_e, edges_e)
    lookup_err_b = lookup_tools.dense_lookup.dense_lookup(err_b, edges_b)
    lookup_err_e = lookup_tools.dense_lookup.dense_lookup(err_e, edges_e)

    # Calculate Et and the Barrel mask
    eta = electrons.eta
    eta_flat = ak.flatten(eta)
    Et_flat = ak.flatten(electrons.pt * (electrons.scEtOverPt + 1))
    eta_barrel_mask = abs(eta_flat) < 1.4442

    # Nominal
    sf_nom_flat = ak.where(
        eta_barrel_mask, 
        lookup_sf_b(Et_flat), 
        lookup_sf_e(Et_flat)
    )
    # Error (Sigma)
    sf_err_flat = ak.where(
        eta_barrel_mask, 
        lookup_err_b(Et_flat), 
        lookup_err_e(Et_flat)
    )

    # Unflatten and attach to electrons
    electrons['trig_eff_ele_nom']  = ak.unflatten(sf_nom_flat, ak.num(eta))
    electrons['trig_eff_ele_up']   = ak.unflatten(sf_nom_flat + sf_err_flat, ak.num(eta))
    electrons['trig_eff_ele_down'] = ak.unflatten(sf_nom_flat - sf_err_flat, ak.num(eta))

    # Fill placeholder muon efficiencies
    electrons['trig_MCeff_mu_nom'] = ak.ones_like(electrons.pt)
    electrons['trig_DATAeff_mu_nom'] = ak.ones_like(electrons.pt)
    electrons['trig_MCeff_mu_up'] = ak.ones_like(electrons.pt)
    electrons['trig_DATAeff_mu_up'] = ak.ones_like(electrons.pt)
    electrons['trig_MCeff_mu_down'] = ak.ones_like(electrons.pt)
    electrons['trig_DATAeff_mu_down'] = ak.ones_like(electrons.pt)
    
    return electrons


def AttachMuonTrigEff(muons, year):
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
    ceval_HLT = correctionlib.CorrectionSet.from_file(ttbarEFT_path(f"data/POG/MUO/{clib_year}/ScaleFactors_Muon_highPt_HLT_{clib_year}_schemaV2.json"))

    # NUM_HLT_DEN_HighPtLooseRelIsoProbes (SF), NUM_HLT_DEN_HighPtLooseRelIsoProbes_MCeff, NUM_HLT_DEN_HighPtLooseRelIsoProbes_DATAeff
    MCeff_nom_flat = ceval_HLT["NUM_HLT_DEN_HighPtLooseRelIsoProbes_MCeff"].evaluate(abseta_flat, pt_flat, "nominal")
    DATAeff_nom_flat = ceval_HLT["NUM_HLT_DEN_HighPtLooseRelIsoProbes_DATAeff"].evaluate(abseta_flat, pt_flat, "nominal")

    MCeff_up_flat = ceval_HLT["NUM_HLT_DEN_HighPtLooseRelIsoProbes_MCeff"].evaluate(abseta_flat, pt_flat, "systup")
    DATAeff_up_flat = ceval_HLT["NUM_HLT_DEN_HighPtLooseRelIsoProbes_DATAeff"].evaluate(abseta_flat, pt_flat, "systup")   

    MCeff_down_flat = ceval_HLT["NUM_HLT_DEN_HighPtLooseRelIsoProbes_MCeff"].evaluate(abseta_flat, pt_flat, "systdown")
    DATAeff_down_flat = ceval_HLT["NUM_HLT_DEN_HighPtLooseRelIsoProbes_DATAeff"].evaluate(abseta_flat, pt_flat, "systdown")    

    muons['trig_MCeff_mu_nom'] = ak.unflatten(MCeff_nom_flat, ak.num(pt))
    muons['trig_DATAeff_mu_nom'] = ak.unflatten(DATAeff_nom_flat, ak.num(pt))
    muons['trig_MCeff_mu_up'] = ak.unflatten(MCeff_up_flat, ak.num(pt))
    muons['trig_DATAeff_mu_up'] = ak.unflatten(DATAeff_up_flat, ak.num(pt))
    muons['trig_MCeff_mu_down'] = ak.unflatten(MCeff_down_flat, ak.num(pt))
    muons['trig_DATAeff_mu_down'] = ak.unflatten(DATAeff_down_flat, ak.num(pt))

    # fill ele trig eff with ones 
    muons['trig_eff_ele_nom'] = ak.ones_like(pt)
    muons['trig_eff_ele_up'] = ak.ones_like(pt)
    muons['trig_eff_ele_down'] = ak.ones_like(pt)


def GetTrigSF(events, lep_cat):

    # Select events based on lepton category
    mask = events[f'is_{lep_cat}']

    leps = ak.pad_none(events.leps_pt_sorted, 2)
    l0 = leps[:,0]
    l1 = leps[:,1]

    def calculate_trigSF_mm(m1, m2, var):
        '''
        m1 : muon1
        m2 : muon2
        var: nom, up, down
        '''
        ed1 = getattr(m1, f"trig_DATAeff_mu_{var}") # eff data muon1
        ed2 = getattr(m2, f"trig_DATAeff_mu_{var}") # eff data muon2
        em1 = getattr(m1, f"trig_MCeff_mu_{var}")   # eff MC muon1
        em2 = getattr(m2, f"trig_MCeff_mu_{var}")   # eff MC muon2

        DATA_eff = 1-(1-ed1)*(1-ed2)
        MC_eff   = 1-(1-em1)*(1-em2)

        return DATA_eff / MC_eff


    def calculate_trigSF_em(l0, l1, var):
        '''
        l0 : ele or mu
        l1 : ele or mu
        var: nom, up, down
        '''        
        DATA_eff = getattr(l0, f"trig_DATAeff_mu_{var}") * getattr(l1, f"trig_DATAeff_mu_{var}")
        MC_eff = getattr(l0, f"trig_MCeff_mu_{var}") * getattr(l1, f"trig_MCeff_mu_{var}")

        return DATA_eff / MC_eff

    # calculate trigger SFs from already saved efficiencies
    if lep_cat == 'ee': 
        calc_nom = l0.trig_eff_ele_nom * l1.trig_eff_ele_nom
        calc_up = l0.trig_eff_ele_up * l1.trig_eff_ele_up
        calc_down = l0.trig_eff_ele_down * l1.trig_eff_ele_down
        
    elif lep_cat == 'mm': 
        calc_nom = calculate_trigSF_mm(l0, l1, "nom")
        calc_up = calculate_trigSF_mm(l0, l1, "up")
        calc_down = calculate_trigSF_mm(l0, l1, "down")

    elif lep_cat == 'em': 
        calc_nom = calculate_trigSF_em(l0, l1, "nom")
        calc_up = calculate_trigSF_em(l0, l1, "up")
        calc_down = calculate_trigSF_em(l0, l1, "down")

    nom = ak.where(mask, calc_nom, 1.0)
    up = ak.where(mask, calc_up, 1.0)
    down = ak.where(mask, calc_down, 1.0)

    return nom, up, down
