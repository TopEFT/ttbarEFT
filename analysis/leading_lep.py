from topcoffea.modules.utils import get_list_of_wc_names
from topcoffea.modules.histEFT import HistEFT
from coffea.nanoevents import NanoAODSchema
from utils.buildLikelihood import full_likelihood
from analysis_tools import genObjectSelection, genEventSelection, isClean
from coffea import processor
from hist import Hist

import awkward as ak
import numpy as np
import torch
import hist

NanoAODSchema.warn_missing_crossrefs = False

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=['ctq8'], 
                  hist_lst = None, dtype=np.float32, 
                  do_errors=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst

        self._dtype = dtype
        self._do_errors = do_errors

        self._likelihood = full_likelihood('/scratch365/cmcgrad2/results/ctq8/linear/-5.0/baseline.yaml')
        
        print("\n\n")
        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)
        print("\n\n")

        # Create the histograms
        self._histo_dict = {
            "genDNNDisc": HistEFT(hist.axis.StrCategory(["genDNNDisc"], name="cat"), 
                                hist.axis.Regular(
                                    start = 0,
                                    stop  = 2400,
                                    bins  = 20, 
                                    name  = "discriminator",
                                    flow  = True
                                ),
                                    wc_names=self._wc_names_lst
                               )
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

    def accumulator(self):
        return self._accumulator
        
    def process(self, events):     

        leps, jets = genObjectSelection(events)
        eventMask  = genEventSelection(leps, jets)
        
        eftCoeffs = events.EFTfitCoefficients[eventMask]

        hout = self._histo_dict

        hout["genDNNDisc"].fill(discriminator=leps.pt[eventMask,0],
                          eft_coeff=eftCoeffs,
                          cat="genDNNDisc"
                         )
        
        return hout

    def postprocess(self, accumulator):
        return accumulator
