from analysis_tools import genEventSelection, genObjectSelection, TensorAccumulator
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

        leps, jets = genObjectSelection(events)
        event_selection_mask = genEventSelection(leps, jets)

        leps  = leps[event_selection_mask]
        jets  = jets[event_selection_mask]
        njets = ak.num(jets)
    
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None
        eft_coeffs = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None


        features = features.concat(torch.from_numpy(np.concatenate([[leps.pt[:,0].to_numpy()], 
                                                                    [leps.pt[:,1].to_numpy()], 
                                                                    [leps.eta[:,0].to_numpy()], 
                                                                    [leps.eta[:,1].to_numpy()], 
                                                                    [leps.phi[:,0].to_numpy()], 
                                                                    [leps.phi[:,1].to_numpy()],  
                                                                    [njets.to_numpy()], 
                                                                    [jets.pt[:,0].to_numpy()], 
                                                                    [jets.pt[:,1].to_numpy()],
                                                                   ]]).T))
    
        fit_coefs = fit_coefs.concat(torch.from_numpy(eft_coeffs))

        return {'features': features, 'fit_coefs': fit_coefs}
    
    def postprocess(self, accumulator):
        return accumulator