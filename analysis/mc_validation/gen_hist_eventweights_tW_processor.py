#!/usr/bin/env python
import awkward as ak
import numpy as np
np.seterr(divide='ignore', invalid='ignore', over='ignore')
import hist
from hist import Hist
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = False

from coffea import processor
from coffea.analysis_tools import PackedSelection

#### TW WC POINTS ####
SM_pt = {"ctGIm": 0.0, "ctGRe":0.0, "cHQ3":0.0, "ctWRe":0.0, "cleQt3Re":0.0,
        "cleQt1Re":0.0, "cQl3":0.0, "cbWRe":0.0, "cHtbRe":0.0}

pt1 = {"ctGIm": -0.4,
        "ctGRe": -0.4,
        "cHQ3" : 1.0,
        "ctWRe" : -1.0,
        "cleQt3Re" : 8.0,
        "cleQt1Re" : 15.0,
        "cQl3"     : 10.0,
        "cbWRe"    : -5.0,
        "cHtbRe"   : 5.0}

pt2 = {"ctGIm": -0.4,
        "ctGRe": -0.4,
        "cHQ3" : 1.0,
        "ctWRe" : -1.0,
        "cleQt3Re" : 5.0,
        "cleQt1Re" : 8.0,
        "cQl3"     : 5.0,
        "cbWRe"    : -3.0,
        "cHtbRe"   : 3.0}

pt3 = {"ctGIm": -0.5,
        "ctGRe": -0.5,
        "cHQ3" : 1.5,
        "ctWRe" : -1.5,
        "cleQt3Re" : 10.0,
        "cleQt1Re" : 18.0,
        "cQl3"     : 12.0,
        "cbWRe"    : -6.0,
        "cHtbRe"   : 6.0}

pt4 = {"ctGIm": -1.5,
        "ctGRe": -1.5,
        "cHQ3" : 2.5,
        "ctWRe" : -2.5,
        "cleQt3Re" : 15.0,
        "cleQt1Re" : 20.0,
        "cQl3"     : 20.0,
        "cbWRe"    : -10.0,
        "cHtbRe"   : 10.0}

new1 = {"ctGIm": -0.2,
        "ctGRe": -0.2,
        "cHQ3" : 1.5,
        "ctWRe" : -1.0,
        "cleQt3Re" : 12.0,
        "cleQt1Re" : 15.0,
        "cQl3"     : 4.0,
        "cbWRe"    : -10.0,
        "cHtbRe"   : 4.0}

new2 = {"ctGIm": -0.4,
        "ctGRe": -0.4,
        "cHQ3" : 1.5,
        "ctWRe" : -1.0,
        "cleQt3Re" : 12.0,
        "cleQt1Re" : 15.0,
        "cQl3"     : 6.0,
        "cbWRe"    : -12.0,
        "cHtbRe"   : 12.0}


# OG stpt
#pt3 = {"ctGIm": -0.4, "ctGRe":-0.4, "cHQ3":1.0, "ctWRe":-1.0, "cleQt3Re":8.0,
#        "cleQt1Re":15.0, "cQl3":10.0, "cbWRe":-5.0, "cHtbRe":5.0}

# SMstpt ("small start point")
#pt1 = {"ctGIm": -0.5, "ctGRe":-0.5, "cHQ3":1.0, "ctWRe":-1.0, "cleQt3Re":1.5,
#        "cleQt1Re":1.5, "cQl3":1.5, "cbWRe":-1.5, "cHtbRe":1.5}

# MED stpt
#pt2 = {"ctGIm": -0.4, "ctGRe":-0.4, "cHQ3":1.0, "ctWRe":-1.0, "cleQt3Re":2.0,
#        "cleQt1Re":4.0, "cQl3":3.0, "cbWRe":-2.0, "cHtbRe":2.0}

# MLG stpt
#pt4 = {"ctGIm": -0.5, "ctGRe":-0.5, "cHQ3":1.5, "ctWRe":-1.5, "cleQt3Re":10.0,
#        "cleQt1Re":10.0, "cQl3":10.0, "cbWRe":-5.0, "cHtbRe":5.0}

# LG stpt
#pt5 = {"ctGIm": -0.5, "ctGRe":-0.5, "cHQ3":2.0, "ctWRe":-2.0, "cleQt3Re":15.0,
#        "cleQt1Re":20.0, "cQl3":20.0, "cbWRe":-10.0, "cHtbRe":10.0}

# Get the lumi for the given year
def get_lumi(year):
    lumi_dict = {
        "2016APV": 19.52,
        "2016": 16.81,
        "2017": 41.48,
        "2018": 59.83
    }
    if year not in lumi_dict.keys():
        raise Exception(f"(ERROR: Unknown year \"{year}\".")
    else:
        return(lumi_dict[year])

# Clean the objects
def is_clean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)

# Create list of wc values in the correct order
def order_wc_values(wcs, ref_pts):
    '''Returns list of wc values in the same order as the list of wc names based on a dictionary
    '''
    wc_names = wcs
    ref_pts = ref_pts

    wc_values = []
    for i in wc_names:
        wc_values.append(ref_pts[i])

    return wc_values

# Calculate event weights from wc values and eft fit coefficients
def calc_event_weights(eft_coeffs, wc_vals):
    '''Returns an array that contains the event weight for each event.
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

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst = None, dtype=np.float32, do_errors=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst

        self._dtype = dtype
        self._do_errors = do_errors

        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)

        # Create the histograms with new scikit hist
        self._histo_dict = {
            "weights_SM"        : Hist(hist.axis.Regular(bins=20, start=-10, stop=10, name="weights_SM", label='event weights at the SM'), storage="weight"),
            "weights_SM_log"    : Hist(hist.axis.Regular(bins=70, start=-6, stop=2, name="weights_SM_log", label='log(event weights at the SM)'), storage="weight"),
            "weights_pt1_log"  : Hist(hist.axis.Regular(bins=70, start=-6, stop=2, name="weights_pt1_log", label='log(event weights at the pt1)'), storage="weight"),
            "weights_pt2_log"  : Hist(hist.axis.Regular(bins=70, start=-6, stop=2, name="weights_pt2_log", label='log(event weights at the pt2)'), storage="weight"),
            "weights_pt3_log"  : Hist(hist.axis.Regular(bins=70, start=-6, stop=2, name="weights_pt3_log", label='log(event weights at the pt3)'), storage="weight"),
            "weights_pt4_log"  : Hist(hist.axis.Regular(bins=70, start=-6, stop=2, name="weights_pt4_log", label='log(event weights at the pt4)'), storage="weight"),
            "weights_new1_log"  : Hist(hist.axis.Regular(bins=70, start=-6, stop=2, name="weights_new1_log", label='log(event weights at new1)'), storage="weight"),
            "weights_new2_log"  : Hist(hist.axis.Regular(bins=70, start=-6, stop=2, name="weights_new2_log", label='log(event weights at new2)'), storage="weight"),

        }

        # Set the list of hists to to fill
        if hist_lst is None:
            self._hist_lst = list(self._histo_dict.keys())
        else:
            for h in hist_lst:
                if h not in self._histo_dict.keys():
                    raise Exception(f"Error: Cannot specify hist \"{h}\", it is not defined in self._histo_dict")
            self._hist_lst = hist_lst

        print("hist_lst: ", self._hist_lst)
    
    @property
    def columns(self):
        return self._columns

    def process(self, events):

        # Dataset parameters
        dataset = events.metadata['dataset']
        hist_axis_name = self._samples[dataset]["histAxisName"]

        year   = self._samples[dataset]['year']
        xsec   = self._samples[dataset]['xsec']
        sow    = self._samples[dataset]['nSumOfWeights']

        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None

        # else: 
        wc_lst_SM = order_wc_values(self._wc_names_lst, SM_pt)
        event_weights_SM = calc_event_weights(eft_coeffs, wc_lst_SM)
            
        wc_lst_pt1 = order_wc_values(self._wc_names_lst, pt1)
        event_weights_pt1 = calc_event_weights(eft_coeffs, wc_lst_pt1)

        wc_lst_pt2 = order_wc_values(self._wc_names_lst, pt2)
        event_weights_pt2 = calc_event_weights(eft_coeffs, wc_lst_pt2)

        wc_lst_pt3 = order_wc_values(self._wc_names_lst, pt3)
        event_weights_pt3 = calc_event_weights(eft_coeffs, wc_lst_pt3)

        wc_lst_pt4 = order_wc_values(self._wc_names_lst, pt4)
        event_weights_pt4 = calc_event_weights(eft_coeffs, wc_lst_pt4)

        wc_lst_new1 = order_wc_values(self._wc_names_lst, new1)
        event_weights_new1 = calc_event_weights(eft_coeffs, wc_lst_new1)

        wc_lst_new2 = order_wc_values(self._wc_names_lst, new2)
        event_weights_new2 = calc_event_weights(eft_coeffs, wc_lst_new2)

        #wc_lst_pt5 = order_wc_values(self._wc_names_lst, pt5)
        #event_weights_pt5 = calc_event_weights(eft_coeffs, wc_lst_pt5)

        counts = np.ones_like(event_weights_SM)

        ######## Fill histos ########
        hout = self._histo_dict
        variables_to_fill = {
            "weights_SM"    : event_weights_SM,
            "weights_SM_log": np.log10(event_weights_SM),
            "weights_pt1_log": np.log10(event_weights_pt1),
            "weights_pt2_log": np.log10(event_weights_pt2),
            "weights_pt3_log": np.log10(event_weights_pt3),
            "weights_pt4_log": np.log10(event_weights_pt4),
            "weights_new1_log": np.log10(event_weights_new1),
            "weights_new2_log": np.log10(event_weights_new2),
            #"weights_pt5_log": np.log10(event_weights_pt5),
        }

        for var_name, var_values in variables_to_fill.items():
            if var_name not in self._hist_lst:
                print(f"Skipping \"{var_name}\", it is not in the list of hists to include")
                continue

            hout[var_name].fill(var_values, weight=counts)

        return hout

    def postprocess(self, accumulator):
        return accumulator

