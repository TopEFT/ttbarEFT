from analysis_tools import event_selection, TensorAccumulator
from coffea.nanoevents import NanoAODSchema
from coffea import processor

import awkward as ak
import numpy as np
import torch

NanoAODSchema.warn_missing_crossrefs = False

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, dtype=torch.float64):
        self._dtype   = dtype
        self._samples = samples
       
    def accumulator(self):
        return self._accumulator
        
    def process(self, events):     

        features  = TensorAccumulator(torch.tensor([]), dtype=self._dtype)
        fit_coefs = TensorAccumulator(torch.tensor([]), dtype=self._dtype)
        
        leps, jets, njets, eft_coeffs = event_selection(events)

        features = features.concat(torch.from_numpy(np.concatenate([[leps.pt[:,0].to_numpy()], 
                                                                    [leps.pt[:,1].to_numpy()], 
                                                                    [leps.eta[:,0].to_numpy()], 
                                                                    [leps.eta[:,1].to_numpy()], 
                                                                    [leps.phi[:,0].to_numpy()], 
                                                                    [leps.phi[:,1].to_numpy()],  
                                                                    [njets.to_numpy()], 
                                                                    [jets.pt[:,0].to_numpy()], 
                                                                    [jets.pt[:,1].to_numpy()],
                                                                    [ak.sum(jets.pt, axis=1).to_numpy()]]).T))
    
        fit_coefs = fit_coefs.concat(torch.from_numpy(eft_coeffs))

        return {'features': features, 'fit_coefs': fit_coefs}
    
    def postprocess(self, accumulator):
        return accumulator