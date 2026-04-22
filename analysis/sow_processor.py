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
import ttbarEFT.modules.corrections as tt_cor 
import topcoffea.modules.corrections as tc_cor
from ttbarEFT.modules.processor_tools import calc_eft_weights


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
            "sow_norm": HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="sow_norm", label="normalized sum of weights for all events"), 
                            wc_names=wc_names_lst, 
                            label="Events"), 
            "nEvents":  HistEFT(
                            proc_axis,
                            hist.axis.Regular(bins=1, start=0, stop=2, name="nEvents", label="number of events"), 
                            wc_names=wc_names_lst, 
                            label="Events"),

            "SumOfWeights":                 HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_ISRUp":           HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_ISRDown":         HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_FSRUp":           HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_FSRDown":         HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),

            "SumOfWeights_renormUp":        HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_renormDown":      HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_factUp":          HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_factDown":        HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_renormfactUp":    HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_renormfactDown":  HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
            "SumOfWeights_toppt":           HistEFT(proc_axis, weights_axis, wc_names=wc_names_lst),
        }

        self._histo_dict['sow_LHEPDFweights'] = hist.Hist(
                proc_axis, 
                weights_axis, 
                hist.axis.Integer(0, 103, name="PDFindex", label="LHEPDFweight Index"),
                storage=hist.storage.Double()
            )


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
        wgts = np.ones_like(events['event'])
        if eft_coeffs is None:
            # Basically any central MC samples
            wgts = events["genWeight"]
        
        norm = xsec/sow

        # attach PS and Qscale weights to events object 
        # tc_cor.AttachPSWeights(events)
        # tt_cor.AttachScaleWeights(events)

        # get arrays of top pt reweights
        LOtoNLO_weights = tt_cor.GetNLO_Weight(events, dataset)
        NLOtoNNLO_weights = tt_cor.GetNNLO_EventWeight(events, dataset)


        ### LHEPdfWeights ###
        if eft_coeffs is not None:
            event_weights_SM = calc_eft_weights(eft_coeffs,np.zeros(len(self._wc_names_lst)))
            pdf_weights = (events.LHEPdfWeight*wgts*event_weights_SM)
        else: 
            pdf_weights = (events.LHEPdfWeight*wgts)
        pdf_index = ak.local_index(pdf_weights, axis=1)

        counts_stacked, index_stacked = ak.broadcast_arrays(counts, pdf_index)


        ####### Fill Histograms #######
        hout = self._histo_dict

        sow_norm_fill_info = {
            "sow_norm"  : counts, 
            "process"   : dataset, 
            "weight"    : wgts*norm, 
            "eft_coeff" : eft_coeffs,
        }
        hout["sow_norm"].fill(**sow_norm_fill_info)

        # Here, weight = counts instead of wgts because this hist is just counting the number
        # of raw events in the file, not effected by the weighting of the event
        nevents_fill_info = {
            "nEvents"   : counts, 
            "process"   : dataset, 
            "weight"    : counts,
            "eft_coeff" : None,
        }
        hout["nEvents"].fill(**nevents_fill_info)

        # Nominal
        hout["SumOfWeights"].fill(process=dataset, SumOfWeights=counts, weight=wgts, eft_coeff=eft_coeffs) #, eft_err_coeff=eft_w2_coeffs)
        
        # Fill ISR/FSR histos
        hout["SumOfWeights_ISRUp"].fill(process=dataset,   SumOfWeights=counts, weight=wgts*events.ISRUp,   eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_ISRDown"].fill(process=dataset, SumOfWeights=counts, weight=wgts*events.ISRDown, eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_FSRUp"].fill(process=dataset,   SumOfWeights=counts, weight=wgts*events.FSRUp,   eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_FSRDown"].fill(process=dataset, SumOfWeights=counts, weight=wgts*events.FSRDown, eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs) 
       
        # Fill renorm/fact histos
        hout["SumOfWeights_renormUp"].fill(process=dataset,       SumOfWeights=counts, weight=wgts*events.renormUp,       eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_renormDown"].fill(process=dataset,     SumOfWeights=counts, weight=wgts*events.renormDown,     eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_factUp"].fill(process=dataset,         SumOfWeights=counts, weight=wgts*events.factUp,         eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_factDown"].fill(process=dataset,       SumOfWeights=counts, weight=wgts*events.factDown,       eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_renormfactUp"].fill(process=dataset,   SumOfWeights=counts, weight=wgts*events.renormfactUp,   eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)
        hout["SumOfWeights_renormfactDown"].fill(process=dataset, SumOfWeights=counts, weight=wgts*events.renormfactDown, eft_coeff=eft_coeffs) # , eft_err_coeff=eft_w2_coeffs)        
        hout["SumOfWeights_toppt"].fill(process=dataset, SumOfWeights=counts, weight=wgts*LOtoNLO_weights*NLOtoNNLO_weights, eft_coeff=eft_coeffs)

        hout['sow_LHEPDFweights'].fill(
            SumOfWeights=ak.flatten(counts_stacked),
            process= dataset,
            PDFindex=ak.flatten(index_stacked),
            weight=ak.flatten(pdf_weights),
        )


        return hout


    def postprocess(self, accumulator):
        return accumulator


