from ttbarEFT.modules.analysis_tools import genEventSelection, genObjectSelection, TensorAccumulator, get_lumi
from coffea.nanoevents import NanoAODSchema
from coffea import processor
from torch import from_numpy, Generator, tensor
from torch.utils.data import random_split, TensorDataset

import awkward as ak
import numpy as np

NanoAODSchema.warn_missing_crossrefs = False

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, dtype):
        self._samples = samples
        self._dtype   = dtype
        
    def accumulator(self):
        return self._accumulator
        
    def process(self, events):  
        dataset = events.metadata['dataset']
        xsec    = self._samples[dataset]['xsec']
        sow     = self._samples[dataset]['nSumOfWeights']
        year    = self._samples[dataset]['year']
        
        factor = xsec/sow*get_lumi(year)*1000

        train_features  = TensorAccumulator(tensor([]), dtype=self._dtype)
        train_fit_coefs = TensorAccumulator(tensor([]), dtype=self._dtype)
        test_features   = TensorAccumulator(tensor([]), dtype=self._dtype)
        test_fit_coefs  = TensorAccumulator(tensor([]), dtype=self._dtype)

        leps, jets = genObjectSelection(events)
        event_selection_mask = genEventSelection(leps, jets)

        leps  = leps[event_selection_mask]
        jets  = jets[event_selection_mask]
        njets = ak.num(jets)
    
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, "EFTfitCoefficients") else None
        eft_coeffs = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None
        eft_coeffs *= factor if eft_coeffs is not None else None


        features = from_numpy(np.concatenate([[leps.pt[:,0].to_numpy()], 
                                              [leps.pt[:,1].to_numpy()], 
                                              [leps.eta[:,0].to_numpy()], 
                                              [leps.eta[:,1].to_numpy()], 
                                              [leps.phi[:,0].to_numpy()], 
                                              [leps.phi[:,1].to_numpy()],  
                                              [njets.to_numpy()], 
                                              [jets.pt[:,0].to_numpy()], 
                                              [jets.pt[:,1].to_numpy()],
                                              [(jets[:,0] + jets[:,1] + leps[:,0] + leps[:,1]).mass.to_numpy()],
                                              [jets[:,0].delta_r(jets[:,1]).to_numpy()],
                                              [(jets[:,0].rho/jets[:,1].rho).to_numpy()],
                                              [abs(jets[:,0].eta - jets[:,1].eta).to_numpy()]
                                             ]).T)
    
        fit_coefs = from_numpy(eft_coeffs)

        train, test, _ = random_split(TensorDataset(features, fit_coefs), [0.10, 0.05, 0.85], generator=Generator().manual_seed(42))

        train_features  = train_features.concat(train[:][0])
        train_fit_coefs = train_fit_coefs.concat(train[:][1])
        test_features   = test_features.concat(test[:][0])
        test_fit_coefs  = test_fit_coefs.concat(test[:][1])

        output = {
            'train_features':  train_features, 
            'train_fit_coefs': train_fit_coefs, 
            'test_features':   test_features, 
            'test_fit_coefs':  test_fit_coefs
        }

        return output
    
    def postprocess(self, accumulator):
        return accumulator