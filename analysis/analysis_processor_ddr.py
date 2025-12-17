#!/usr/bin/env python
import copy
import coffea
import numpy as np
import awkward as ak
import json
import hist
import yaml

# from mt2 import mt2

from coffea import processor
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods import vector
from coffea.lumi_tools import LumiMask

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

from topcoffea.modules.paths import topcoffea_path
from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth
import topcoffea.modules.corrections as tc_cor
import topcoffea.modules.event_selection as tc_es

from ttbarEFT.modules.paths import ttbarEFT_path
from ttbarEFT.modules.analysis_tools import make_mt2, get_lumi
import ttbarEFT.modules.object_selection as tt_os
import ttbarEFT.modules.event_selection as tt_es
from ttbarEFT.modules.corrections import AttachElectronSF, AttachMuonSF, AttachMuonTrigSF, ApplyMuonPtCorr, ApplyJetVetoMaps

from topcoffea.modules.get_param_from_jsons import GetParam

get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_tt_param = GetParam(ttbarEFT_path("params/params.json"))

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, lep_cat, wc_names_lst=[], hist_lst=None, dtype=np.float32):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 
        self._lep_cat = lep_cat

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)
        # appl_axis = hist.axis.StrCategory([], name="appl", label=r"AR/SR", growth=True)

        # fill histograms using info from axes.json
        with open(ttbarEFT_path("params/axes.json"), 'r') as axes_file:
            axes_info = json.load(axes_file)

        histograms = {}
        for name, info in axes_info['CR_axes'].items():
            if 'variable' in info: 
                dense_axis = hist.axis.Variable(info['variable'], name=name, label=info['label'])
                sumw2_axis = hist.axis.Variable(info['variable'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')
            else:
                dense_axis = hist.axis.Regular(*info['regular'], name=name, label=info['label'])
                sumw2_axis = hist.axis.Regular(*info['regular'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')

            histograms[name] = HistEFT(
                proc_axis, 
                # chan_axis,
                # syst_axis,
                dense_axis,
                wc_names = wc_names_lst, 
                label=r'Events',
            )

            # histograms[name+'_sumw2'] = HistEFT(
            #     proc_axis, 
            #     # chan_axis,
            #     # syst_axis,
            #     sumw2_axis,
            #     wc_names = wc_names_lst, 
            #     label=r'Events',
            # )

        self._accumulator = histograms

        # set the list of hists to fill
        if hist_lst is None:
            self._hist_lst = list(self._accumulator.keys()) #fill all hists if not specified
        else:
            for hist_to_include in hist_lst:
                if hist_to_include not in self._accumulator.keys():
                    raise Exception(f"Error: Cannot specify hist \'{hist_to_include}\', it is not defined in the processor.")
            self._hist_lst = hist_lst 

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    def process(self, events):

        # Dataset parameters
        dataset         = events.metadata['dataset']
        isData          = self._samples[dataset]['isData']
        histAxisName    = self._samples[dataset]['histAxisName']
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']

        print(f"\n\n histAxisName: {histAxisName}")
        print(f"dataset: {dataset}")
        print(f"year: {year}")
        print(f"xsec: {xsec} \n\n")

        isEFT = hasattr(events, 'EFTfitCoefficients')    
        assert not (isEFT and isData), f"isEFT and isData cannot both be True. Check input samples."

        lep_cat = self._lep_cat

        datasets = ["Muon", "SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG"]
        for d in datasets:
            if dataset.startswith(d):
                dataset = dataset.split('_')[0]


        ######### EFT coefficients ##########
        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, 'EFTfitCoefficients') else None
        if eft_coeffs is not None:
            # Check to see if the ordering of WCs for this sample matches what want
            if self._samples[dataset]['WCnames'] != self._wc_names_lst:
                eft_coeffs = efth.remap_coeffs(self._samples[dataset]['WCnames'], self._wc_names_lst, eft_coeffs)

        # Initialize the out object
        hout = self.accumulator


        ######### Create Lepton Categories ##########
        cat_dict = None
        with open(ttbarEFT_path("params/channels.yaml"), "r") as f:
            cat_dict=yaml.safe_load(f)

        CR_cat_dict = cat_dict['CR_CHANNELS']
        SR_cat_dict = cat_dict['SR_CHANNELS']


        ######### Initialize Objects #########
        met  = events.MET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet 

        leptonSelection = tt_os.Run2LeptonSelection()

        # An array of length events that is just 1 for each event
        events['nom'] = ak.ones_like(events.MET.pt)

        ######### Electron Selection ##########
        ele['isGoodElec']=leptonSelection.is_sel_ele(ele)
        ele_good = ele[ele.isGoodElec]
        # ele_good = ele

        # test_pt = ele.pt


        ######### Muon Selection ##########
        mu['pt'] = ApplyMuonPtCorr(mu, year, isData)
        mu['isGoodMuon']=leptonSelection.is_sel_muon(mu)
        mu_good = mu[mu.isGoodMuon]


        ######### Lepton Selection ##########
        # add lepton scale factors 
        if not isData: 
            AttachElectronSF(ele_good, year)    
            AttachMuonSF(mu_good, year)     

        leps = ak.concatenate([ele_good, mu_good], axis=1)
        leps_sorted = leps[ak.argsort(leps.pt, axis=-1,ascending=False)] 

        events['leps_pt_sorted'] = leps_sorted


        ######### Jet Selections #########
        jets['isPres'] = tt_os.is_pres_jet(jets) 
        goodJets = jets[jets.isPres] 

        njets = ak.num(goodJets)
        ht = ak.sum(goodJets.pt,axis=-1)
        j0 = goodJets[ak.argmax(goodJets.pt,axis=-1,keepdims=True)]

        # Medium DeepJet WP
        medium_tag = "btag_wp_medium_" + year.replace("201", "UL1")
        btagwpm = get_tc_param(medium_tag)
        isBtagJetsMedium = (goodJets.btagDeepFlavB > btagwpm)
        isNotBtagJetsMedium = np.invert(isBtagJetsMedium)
        nbtagsm = ak.num(goodJets[isBtagJetsMedium])

        # JetVetoMaps applied
        veto_map_array = ApplyJetVetoMaps(goodJets, year)   
        veto_map_mask = (veto_map_array == 0)


        ######### Systematics #########

        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)

        if not isData: 
            # If this is not an eft sample, get the genWeight
            if eft_coeffs is None: 
                genw = events['genWeight']
            else: 
                genw = np.ones_like(events['event'])

            lumi = 1000.0*get_lumi(year)
            norm = genw*(xsec/sow)*lumi
            weights_obj_base.add('norm', norm)

            tt_es.addLepSFs(events, ele_good, mu_good) 
            # AttachPSWeights(events)
            # AttachScaleWeights(events) 
            # AttachPdfWeights(events)

            #weights_obj_base.add('ISR')...
            #weights_obj_base.add('FSR')...
            #weights_obj_base.add('renorm')...
            #weights_obj_base.add('fact')...
            #weights_obj_base.add('Prefiring')...
            #weights_obj_base.add('PU')...

            # TODO: Are these the correct way to apply these? Might need to be applied based on lepton category
            # weights_obj_base.add('lepSF_ele', events.SF_2l_ele, copy.deepcopy(events.SF_2l_ele_up), copy.deepcopy(events.SF_2l_ele_down))
            # weights_obj_base.add('lepSF_muon', events.SF_2l_muon, copy.deepcopy(events.SF_2l_muon_up), copy.deepcopy(events.SF_2l_muon_down)) 


        ######### Add variables to EVENTS #########
        events['njets'] = ak.num(jets)
        tt_es.addLepCatMasks(events) 

        ######### Create objects for dense axes ##########
        leps_sorted = ak.pad_none(leps_sorted, 2)
        l0 = leps_sorted[:,0]
        l1 = leps_sorted[:,1]

        ptll = (l0+l1).pt
        mll = (l0+l1).mass


        ######### Selection Masks #########

        if isData:
            # Lumi Mask for Data    
            golden_json_path = {
                "2016": topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
                "2016APV": topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
                "2017": topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
                "2018": topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
            }
            lumi_mask = LumiMask(golden_json_path[year])(events.run,events.luminosityBlock)

        # Charge masks
        chargel0_p = ak.fill_none(((l0.charge)>0),False)
        chargel0_m = ak.fill_none(((l0.charge)<0),False)
        charge2l_os = ak.fill_none(((l0.charge+l1.charge)==0),False)
        charge2l_ss = ak.fill_none(((l0.charge+l1.charge)!=0),False)


        ######### Selection Masks #########
        pass_trg = tt_es.trg_pass_no_overlap(events, isData, dataset, str(year), tt_es.triggers_dict, tt_es.exclude_triggers_dict, lep_cat)
        # pass_trg = tc_es.trg_pass_no_overlap(events, isData, dataset, str(year), tt_es.triggers_dict, tt_es.exclude_triggers_dict)
        # at_least_two_leps = ak.fill_none(nleps>=2, False)

        # Charge masks
        chargel0_p = ak.fill_none(((l0.charge)>0),False)
        chargel0_m = ak.fill_none(((l0.charge)<0),False)
        charge2l_os = ak.fill_none(((l0.charge+l1.charge)==0),False)
        charge2l_ss = ak.fill_none(((l0.charge+l1.charge)!=0),False)


        ######### Store boolean masks with PackedSelection ##########
        selections = PackedSelection(dtype='uint64')
        
        if isData: 
            selections.add('is_good_lumi', lumi_mask)

        selections.add('pass_trg', pass_trg)
        selections.add('jetvetomap', veto_map_mask)

        # selections.add('2l', at_least_two_leps) 
        selections.add('2los', charge2l_os)

        selections.add("ee",  (events.is_ee & pass_trg))
        selections.add("em",  (events.is_em & pass_trg))
        selections.add("mm",  (events.is_mm))

        selections.add('bmask_exactly0med', (nbtagsm==0))
        selections.add('bmask_exactly1med', (nbtagsm==1))
        selections.add('bmask_exactly2med', (nbtagsm==2))
        selections.add('bmask_atleast2med', (nbtagsm>=2))

        selections.add("exactly_0j", (njets==0))
        selections.add("exactly_1j", (njets==1))
        selections.add("exactly_2j", (njets==2))

        selections.add("atleast_1j", (njets>=1))

        ######### Fill dense axes variables ##########

        dense_axis_variables = {}
        dense_axis_variables['njets'] = njets
        dense_axis_variables['nbjets'] = nbtagsm
        dense_axis_variables['mll'] = mll
        dense_axis_variables['ptll'] = ptll
        dense_axis_variables['l0pt'] = l0.pt
        dense_axis_variables['l0eta'] = l0.eta 
        dense_axis_variables['l0phi'] = l0.phi 
        dense_axis_variables['l1pt'] = l1.pt
        dense_axis_variables['l1eta'] = l1.eta 
        dense_axis_variables['l1phi'] = l1.phi
        dense_axis_variables['j0pt'] = ak.flatten(j0.pt) 
        dense_axis_variables['j0eta'] = ak.flatten(j0.eta)
        dense_axis_variables['j0phi'] = ak.flatten(j0.phi)
        dense_axis_variables['MET'] = met.pt

        jet_variables = ['j0pt', 'j0eta', 'j0phi']


        ########## Fill the histograms ##########
        # Loop through systematics and fill histograms
        weights_object = copy.deepcopy(weights_obj_base)
        weight = weights_object.weight(None) 
        # if eft_coeffs is None:
            # weight= events["genWeight"]
        # else:
        # weight = np.ones_like(events['event'])

        for jet_cat in CR_cat_dict[lep_cat]['jet_list']: 
            # masks that are applied to all categories
            cuts_list = ['jetvetomap', '2los', 'bmask_exactly0med']

            if isData:
                cuts_list.append('is_good_lumi')
            
            cuts_list.append(lep_cat)
            cuts_list.append(jet_cat)


            event_selection_mask = selections.all(*(cuts_list))
            eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

            for dense_axis_name, dense_axis_vals in dense_axis_variables.items():
                # if the category requires zero jets, don't fill jet histograms
                if (jet_cat == 'exactly_0j') and (dense_axis_name in jet_variables):
                    print(f"Skipping '{dense_axis_name}' in category '{lep_cat}_{jet_cat}'. Jet histograms are not filled for categories that don't require a jet")
                    continue

                # if dense_axis_name not in self._hist_lst:
                #     print(f"Skipping \"{dense_axis_name}\", it is not in the list of hists to include")
                #     continue                         

                # Fill the histos
                axes_fill_info_dict = {
                    dense_axis_name : dense_axis_vals[event_selection_mask],
                    "process"       : histAxisName,
                    "weight"        : weight[event_selection_mask],
                    "eft_coeff"     : eft_coeffs_cut,
                }

                hout[dense_axis_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator


