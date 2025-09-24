from awkward import num
from coffea.nanoevents.NanoAODSchema import warn_missing_crossrefs
from coffea import processor
from numpy import concatenate
from os import makedirs
from random import randint
from time import time
from torch import from_numpy, Generator, tensor, save
from torch.utils.data import random_split, TensorDataset
from ttbarEFT.modules.analysis_tools import genEventSelection, genObjectSelection, get_lumi

warn_missing_crossrefs = False

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, dtype, wcs, outname):
        self._samples = samples
        self._dtype   = dtype
        self._wcs     = wcs
        self._outname = outname 

        makedirs(f'{outname}/train', mode=0o755, exist_ok=True)
        makedirs(f'{outname}/test',  mode=0o755, exist_ok=True)
        
    def accumulator(self):
        return self._accumulator
        
    def process(self, events):
        dataset = events.metadata['dataset']
        xsec    = self._samples[dataset]['xsec']
        sow     = self._samples[dataset]['nSumOfWeights']
        year    = self._samples[dataset]['year']
        
        factor = xsec/sow*get_lumi(year)*1000

        leps, jets = genObjectSelection(events)
        event_selection_mask = genEventSelection(leps, jets)

        leps  = leps[event_selection_mask]
        jets  = jets[event_selection_mask]
        pMask = leps.pdgId > 0
        nMask = leps.pdgId < 0
        njets = num(jets)
    
        eft_coeffs = events['EFTfitCoefficients'].to_numpy() if hasattr(events, "EFTfitCoefficients") else None
        eft_coeffs = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None
        eft_coeffs *= factor if eft_coeffs is not None else None

        interference = [int(0.5*x**2+0.5*x) for x in range(len(self._wcs) + 1)]
        quadratic    = [int(0.5*x**2+1.5*x) for x in range(len(self._wcs) + 1)]
        
        for i in range(len(interference)):
            if (i == 0) & (len(eft_coeffs[eft_coeffs[:,0] < 0] > 0)): 
                eft_coeffs[eft_coeffs[:,0] < 0, interference] = 0
            elif len(eft_coeffs[eft_coeffs[:,quadratic[i]] < 0]) > 0: 
                eft_coeffs[eft_coeffs[:,quadratic[i]] < 0, interference[i]:quadratic[i]+1] = 0

        features = from_numpy(concatenate([[leps.pt[pMask][:,0].to_numpy()], 
                                           [leps.pt[nMask][:,0].to_numpy()],
                                           [leps.eta[pMask][:,0].to_numpy()],
                                           [leps.eta[nMask][:,0].to_numpy()],
                                           [leps.phi[pMask][:,0].to_numpy()],
                                           [leps.phi[nMask][:,0].to_numpy()],
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
        _ = None

        save(TensorDataset(train[:][0], train[:][1]), f'{self._outname}/train/{int(time())}{randint(10,99)}.p')
        save(TensorDataset(test[:][0],  test[:][1]),  f'{self._outname}/test/{int(time())}{randint(10,99)}.p')

    
    def postprocess(self, accumulator):
        return accumulator