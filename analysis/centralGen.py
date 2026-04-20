from awkward import num
from coffea.nanoevents import NanoAODSchema 
from coffea import processor
from numpy import array, concatenate, float64, float32, sqrt
from os import makedirs
from random import randint
from time import time
from torch import from_numpy, Generator, tensor, save, where
from torch.utils.data import random_split, TensorDataset
from ttbarEFT.modules.analysis_tools import genEventSelection, genObjectSelection, get_lumi, TensorAccumulator

NanoAODSchema.warn_missing_crossrefs = False

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, dtype, wcs, outname):
        self._samples = samples
        self._dtype   = dtype
        self._wcs     = wcs
        self._outname = outname 

        makedirs(f'{outname}/train',       mode=0o755, exist_ok=True)
        makedirs(f'{outname}/test',        mode=0o755, exist_ok=True)
        
    def accumulator(self):
        return self._accumulator
        
    def process(self, events):
        dataset = events.metadata['dataset']
        xsec    = self._samples[dataset]['xsec']
        sow     = self._samples[dataset]['nSumOfWeights']
        year    = self._samples[dataset]['year']
        nEvents = self._samples[dataset]['nEvents']

        out = TensorAccumulator(tensor([]), dtype=self._dtype)
        
        factor = xsec*get_lumi(year)*1000/sow

        leps, jets = genObjectSelection(events)
        event_selection_mask = genEventSelection(leps, jets)

        leps  = leps[event_selection_mask]
        jets  = jets[event_selection_mask]

        eft_coeffs  = events['EFTfitCoefficients'].to_numpy().astype(self._dtype)
        eft_coeffs  = eft_coeffs[event_selection_mask] 
        eft_coeffs *= factor 
        
        pMask = leps.pdgId < 0
        nMask = leps.pdgId > 0
        dilep = leps.sum()
        dijet = jets[:,0] + jets[:,1]

        interference = [int(0.5*x**2+0.5*x) for x in range(len(self._wcs) + 1)][1:]
        quadratic    = [int(0.5*x**2+1.5*x) for x in range(len(self._wcs) + 1)][1:]

        negSm = eft_coeffs[:,0] <= 0
        if len(negSm) > 0:  # if sm is 0 or negative, interference and sm is 0
            eft_coeffs[negSm, 0] = 0
            for i in range(len(interference)):
                eft_coeffs[negSm, interference[i]] = 0
        
        for i in range(len(self._wcs)):
            negQuad = eft_coeffs[:,quadratic[i]] < 0
            if len(eft_coeffs[negQuad]) > 0:  # if quad is negative , set bsm terms to 0
                eft_coeffs[negQuad, interference[i]:quadratic[i] + 1] = 0

        features = from_numpy(concatenate([[leps.pt[:,0].to_numpy()], #leading lepton
                                           [leps.eta[:,0].to_numpy()],
                                           [leps.phi[:,0].to_numpy()],
                                           [leps.pt[:,1].to_numpy()], #sub-leading lepton
                                           [leps.eta[:,1].to_numpy()],
                                           [leps.phi[:,1].to_numpy()],
                                           [leps.pt[nMask][:,0].to_numpy()], #leading negative lepton0
                                           [leps.eta[nMask][:,0].to_numpy()],
                                           [leps.phi[nMask][:,0].to_numpy()],
                                           [leps.pt[pMask][:,0].to_numpy()], #leading positive lepton 
                                           [leps.eta[pMask][:,0].to_numpy()],
                                           [leps.phi[pMask][:,0].to_numpy()],
                                           [dilep.pt.to_numpy()], #dilep system
                                           [dilep.mass.to_numpy()],
                                           [leps[nMask][:,0].delta_phi(leps[pMask][:,0]).to_numpy()],       # ln_lp
                                           [abs(leps[nMask][:,0].eta - leps[pMask][:,0].eta).to_numpy()],   # ln_lp
                                           [jets.pt[:,0].to_numpy()], #leading jet
                                           [jets.eta[:,0].to_numpy()],
                                           [jets.phi[:,0].to_numpy()],
                                           [jets.pt[:,1].to_numpy()], #sub-leading jet
                                           [jets.eta[:,1].to_numpy()],
                                           [jets.phi[:,1].to_numpy()],
                                           [dijet.pt.to_numpy()], #dijet system
                                           [dijet.mass.to_numpy()],
                                           [jets[:,0].delta_phi(jets[:,1]).to_numpy()], #                  j1_j2
                                           [abs(jets[:,0].eta - jets[:,1].eta).to_numpy()], #              j1_j2
                                           [(dilep + dijet).mass.to_numpy()], # psuedo mtt
                                           [num(jets).to_numpy()],  #                                      nJets
                                           [concatenate([[dijet.pt.to_numpy()], #lj0pt
                                                         [dilep.pt.to_numpy()], 
                                                         [(leps[:,0] + jets[:,0]).pt.to_numpy()]], 
                                                        axis=0).max(0)],
                                          ]).astype(self._dtype).T)
        
        massMask = where(~features[:,13].isnan())
        if massMask[0].sum() > 0:
            print(f'{massMask[0].sum():,} events cut from nan mass values!')
            
        fit_coefs = from_numpy(eft_coeffs)

        train, test = random_split(TensorDataset(features[massMask], fit_coefs[massMask]), [0.75, 0.25], generator=Generator().manual_seed(42))

        save(TensorDataset(train[:][0], train[:][1]), f'{self._outname}/train/{int(time())}{randint(1000000,9999999)}.p')
        save(TensorDataset(test[:][0],  test[:][1]),  f'{self._outname}/test/{int(time())}{randint(1000000,9999999)}.p')

        output = {'out':  out}

        return output

    def postprocess(self, accumulator):
        return accumulator