import numpy as np
import yaml 

from ttbarEFT.modules.paths import ttbarEFT_path
import ttbarEFT.modules.corrections as tt_cor 


def calc_eft_weights(eft_coeffs, wc_vals):
    '''
    Returns an array that contains the event weight for each event.
    eft_coeffs: Array of eft fit coefficients for each event
    wc_vals: wilson coefficient values desired for the event weight calculation, listed in the same order as the wc_lst
             such that the multiplication with eft_coeffs is correct
             The correct ordering can be achieved with the order_wc_values function
    '''
    event_weight = np.empty_like(eft_coeffs)

    wcs = np.hstack((np.ones(1),wc_vals))
    wc_cross_terms = []
    index = 0
    for j in range(len(wcs)):
        for k in range (j+1):
            term = wcs[j]*wcs[k]
            wc_cross_terms.append(term)
    event_weight = np.sum(np.multiply(wc_cross_terms, eft_coeffs), axis=1)

    return event_weight



def get_syst_lists(year, isData, syst_list=[], run_era=None):

    with open(ttbarEFT_path("params/syst_names.yaml"), "r") as f:
        syst_names_yaml=yaml.safe_load(f)

    wgt_correction_bases = syst_names_yaml['wgt_correction_bases']
    btag_var = syst_names_yaml[f"btag_var_{year}"]

    obj_correction_syst_lst = tt_cor.get_supported_jet_systematics(year, isData=isData, era=run_era)
    kinematic_variations = ['nominal']
    event_weight_variations = []

    if (syst_list) and (not isData):                                    # if doing systematics, loop over corrections for only MC
        if 'onlyJEC' in syst_list:                                      # if onlyJEC, just fill the kinematic variations
            kinematic_variations.extend(obj_correction_syst_lst)
            event_weight_variations = []

        elif 'onlyEventWeights' in syst_list:
            for w in wgt_correction_bases:                              # if all event weight systematics, loop through all bases from yaml 
                if w == 'btagSF': 
                    for v in btag_var:
                        event_weight_variations.extend([f"{v}Up"])
                        event_weight_variations.extend([f"{v}Down"])
                else: 
                    event_weight_variations.extend([f"{w}Up"])
                    event_weight_variations.extend([f"{w}Down"])

        elif 'all' in syst_list: 
            kinematic_variations.extend(obj_correction_syst_lst)        # if all systematics, include all object corrections 

            for w in wgt_correction_bases:                              # if all systematics, loop through all bases from yaml 
                if w == 'btagSF': 
                    for v in btag_var:
                        event_weight_variations.extend([f"{v}Up"])
                        event_weight_variations.extend([f"{v}Down"])
                else: 
                    event_weight_variations.extend([f"{w}Up"])
                    event_weight_variations.extend([f"{w}Down"])
        else:                                                           # else, loop through just syst variations in provided list
            for var in syst_list: 
                if var in wgt_correction_bases: 
                    if var == 'btagSF':
                        for v in btag_var:
                            event_weight_variations.extend([f"{v}Up"])
                            event_weight_variations.extend([f"{v}Down"])
                    else:
                        event_weight_variations.extend([f"{var}Up"])
                        event_weight_variations.extend([f"{var}Down"])

                if var in obj_correction_syst_lst: 
                    kinematic_variations.extend([var])                                     

        if 'noJEC' in syst_list:                                        # set syst_var_list to empty if noJEC specified
            kinematic_variations = ['nominal']


    return event_weight_variations, kinematic_variations

