'''
This processor processes a root file and returns a 1 bin histogram that just holds one entry per 
event weighted by the EFTFitCoeff. No cuts or event selections are made.
Taking this histogram, do histEFT.as_hist({}) to get the sow of all events at the SM for the sample. 
This is equivalent to summing EFTFitCoeff[0] for all events for an EFT sample. 
'''

#!/usr/bin/env python
import numpy as np
import awkward as ak

np.seterr(divide='ignore', invalid='ignore', over='ignore')
from coffea import processor

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = False

import hist
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth


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


# Main analysis processor
class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst = None, dtype=np.float32, do_errors=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst

        self._dtype = dtype
        self._do_errors = do_errors

        print("\n\n")
        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)
        print("\n\n")

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        weights_axis = hist.axis.Regular(bins=1, start=0, stop=2, name="SumOfWeights", label="SumOfWeights")

        self._histo_dict = {
            "nEvents":      HistEFT(
                                proc_axis,
                                hist.axis.Regular(bins=1, start=0, stop=2, name="nEvents", label="number of events"), 
                                wc_names=wc_names_lst, 
                                label="Events"),
            "SumOfWeights": HistEFT(
                                proc_axis, 
                                hist.axis.Regular(bins=1, start=0, stop=2, name="SumOfWeights", label="SumOfWeights"), 
                                wc_names=wc_names_lst,
                                label="Events"),
            "SumOfWeights_SM": HistEFT(
                                proc_axis, 
                                hist.axis.Regular(bins=1, start=0, stop=2, name="SumOfWeights_SM", label="SumOfWeights_SM"), 
                                wc_names=wc_names_lst,
                                label="Events"),
            "sumw2":        HistEFT(
                                proc_axis, 
                                hist.axis.Regular(bins=1, start=0, stop=2, name="sumw2", label="sumw2"),
                                wc_names = wc_names_lst,
                                label="Events"),
            "weights_SM":   HistEFT(
                                proc_axis, 
                                hist.axis.Regular(bins=40, start=0, stop=20, name="weights_SM", label='event weights at SM'),
                                wc_names = wc_names_lst,
                                label="Events"),
            "weights_SM_log": HistEFT(
                                proc_axis, 
                                hist.axis.Regular(bins=120, start=-4, stop=2, name="weights_SM_log", label='log(event weights at the SM)'),
                                wc_names = wc_names_lst,
                                label="Events"),

        }

    @property
    def columns(self):
        return self._columns

    def process(self, events):

        # Dataset parameters
        dataset = events.metadata['dataset']
        isData  = self._samples[dataset]["isData"]
        if isData: raise Exception("Why are you running this over data?")

        hist_axis_name = self._samples[dataset]["histAxisName"]
        sow = self._samples[dataset]['nSumOfWeights']
        xsec = self._samples[dataset]['xsec']

        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None
        # eft_w2_coeffs = efth.calc_w2_coeffs(eft_coeffs,self._dtype) if (self._do_errors and eft_coeffs is not None) else None

        # Get nominal wgt
        counts = np.ones_like(events['event'])
        weight = np.ones_like(events['event'])
        if eft_coeffs is None:
            # Basically any central MC samples
            weight = events["genWeight"]
        

        if eft_coeffs is not None:
            event_weights = calc_eft_weights(eft_coeffs,np.zeros(len(self._wc_names_lst)))*weight
            sumw2 = np.square(event_weights)
        else: 
            event_weights = weight
            sumw2 = np.square(weight)


        ####### Fill Histogram #######
        hout = self._histo_dict

        hout['nEvents'].fill(process=dataset, nEvents=counts, weight=counts, eft_coeff=None)
        hout['SumOfWeights'].fill(process=dataset, SumOfWeights=counts, weight=weight, eft_coeff=eft_coeffs)
        hout['SumOfWeights_SM'].fill(process=dataset, SumOfWeights_SM=counts, weight=event_weights, eft_coeff=None)
        hout['sumw2'].fill(process=dataset, sumw2=counts, weight=sumw2, eft_coff=None)
        hout['weights_SM'].fill(process=dataset, weights_SM=event_weights, weight=counts, eft_coeff=None)
        hout['weights_SM_log'].fill(process=dataset, weights_SM_log=np.log10(event_weights), weight=counts, eft_coeff=None)

        return hout


    def postprocess(self, accumulator):
        return accumulator


