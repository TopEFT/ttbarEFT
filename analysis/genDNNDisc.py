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
                                    stop  = 5,
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
        
        leps      = leps[eventMask]
        jets      = jets[eventMask]
        eftCoeffs = events.EFTfitCoefficients[eventMask]
        
        leps, jets = genObjectSelection(events)
        eventMask = genEventSelection(leps, jets)

        #build dict of arrays
        allret={}
        
        for i in range(2):
            allret[f'Lep{i+1}_pt']  = leps.pt[eventMask,i]
            allret[f'Lep{i+1}_eta'] = leps.eta[eventMask,i]
            allret[f'Lep{i+1}_phi'] = leps.phi[eventMask,i]

        for i in range(2):
            allret[f'jet{i+1}_pt']  = jets.pt[eventMask,i]    
        allret['nJet30'] = ak.num(jets[eventMask])
        
        ordered = ['Lep1_pt', 'Lep2_pt', 'Lep1_eta','Lep2_eta', 'Lep1_phi', 'Lep2_phi', 'nJet30', 'jet1_pt', 'jet2_pt']
        allret  = {k: allret[k] for k in ordered}

        kin = torch.tensor([])

        rwt_point = 'ctq8=1.0'

        rwt_point = rwt_point.replace("=","_").replace(":","_").split('_')
        rwt_point = dict(map(lambda i: (rwt_point[i], float(rwt_point[i+1])), range(len(rwt_point)-1)[::2]))

        hout = self._histo_dict
        
        for var_name, var_values in allret.items():
            temp = torch.from_numpy(ak.to_numpy(var_values)).unsqueeze(1)
            kin  = torch.cat((kin,temp), 1)
        
        disc = self._likelihood(kin, rwt_point)
        
        nanMask = (~torch.isnan(disc)).detach().numpy()
        disc = disc[nanMask].detach().numpy()

        hout["genDNNDisc"].fill(discriminator=disc,
                          eft_coeff=eftCoeffs[nanMask],
                          cat="genDNNDisc"
                         )
        
        return hout

    def postprocess(self, accumulator):
        return accumulator
