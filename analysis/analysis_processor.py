#!/usr/bin/env python
import copy
import coffea
import numpy as np
import awkward as ak
import json
import hist

from mt2 import mt2

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
# from ttbarEFT.modules.axes import info as axes_info
import ttbarEFT.modules.object_selection as tt_os
import ttbarEFT.modules.event_selection as tt_es
from ttbarEFT.modules.corrections import ApplyJetVetoMaps

from topcoffea.modules.get_param_from_jsons import GetParam

get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_tt_param = GetParam(ttbarEFT_path("params/params.json"))

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst=None, ecut_threshold=None, do_errors=False, do_systs=False, split_by_lepton_flavor=False, skip_signal_regions=False, skip_control_regions=False, muonSyst='nominal', dtype=np.float32, rebin=False, offZ_split=False, tau_h_analysis=False, fwd_analysis=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst
        self._dtype = dtype 

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        # syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)
        # appl_axis = hist.axis.StrCategory([], name="appl", label=r"AR/SR", growth=True)

        # Set the booleans
        self._do_errors = do_errors # Whether to calculate and store the w**2 coefficients
        self._do_systematics = do_systs # Whether to process systematic samples
        self._split_by_lepton_flavor = split_by_lepton_flavor # Whether to keep track of lepton flavors individually
        self._skip_signal_regions = skip_signal_regions # Whether to skip the SR categories
        self._skip_control_regions = skip_control_regions # Whether to skip the CR categories

        # fill histograms using info from axes.py 
        # histograms = {}
        # for name, info in axes_info.items():
        #     if 'variable' in info: 
        #         dense_axis = hist.axis.Variable(info['variable'], name=name, label=info['label'])
        #         sumw2_axis = hist.axis.Variable(info['variable'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')
        #     else:
        #         dense_axis = hist.axis.Regular(*info['regular'], name=name, label=info['label'])
        #         sum2w_axis = hist.axis.Regular(*info['regular'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')

        #     histograms[name] = HistEFT(
        #         proc_axis, 
        #         syst_axis,
        #         dense_axis,
        #         wc_names = wc_names_lst, 
        #         label=r'Events',
        #     )

        #     histograms[name+'_sumw2'] = HistEFT(
        #         proc_axis, 
        #         syst_axis,
        #         sum2w_axis,
        #         wc_names = wc_names_lst, 
        #         label=r'Events',
        #     )

        with open(ttbarEFT_path("params/axes.json"), 'r') as axes_file:
            axes_info = json.load(axes_file)

        histograms = {}
        for name, info in axes_info['CR_axes'].items():
            if 'variable' in info: 
                dense_axis = hist.axis.Variable(info['variable'], name=name, label=info['label'])
                sumw2_axis = hist.axis.Variable(info['variable'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')
            else:
                dense_axis = hist.axis.Regular(*info['regular'], name=name, label=info['label'])
                sum2w_axis = hist.axis.Regular(*info['regular'], name=name+'_sumw2', label=info['label'] + ' sum of w^2')

            histograms[name] = HistEFT(
                proc_axis, 
                chan_axis,
                # syst_axis,
                dense_axis,
                wc_names = wc_names_lst, 
                label=r'Events',
            )

            # histograms[name+'_sumw2'] = HistEFT(
            #     proc_axis, 
            #     syst_axis,
            #     sum2w_axis,
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

        # print out basic info before running over filesrun2leptonselection
        print("\n\n")
        print("self._samples", self._samples)
        print("self._wc_names_lst", self._wc_names_lst)
        print("hist_lst: ", self._hist_lst)
        print("\n\n")

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns


    def process(self, events):

        # Dataset parameters
        dataset         = events.metadata['dataset']
        # isEFT             = self._samples[dataset]["WCnames"] != []
        isData          = self._samples[dataset]['isData']
        histAxisName  = self._samples[dataset]['histAxisName']
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']

        isEFT = hasattr(events, 'EFTfitCoefficients')    
        assert not (isEFT and isData), f"isEFT and isData cannot both be True. Check input samples."


        ######### EFT coefficients ##########
        # Extract the EFT quadratic coefficients and optionally use them to calculate the coefficients on the w**2 quartic function
        # eft_coeffs is never Jagged so convert immediately to numpy for ease of use.
        eft_coeffs = ak.to_numpy(events['EFTfitCoefficients']) if hasattr(events, 'EFTfitCoefficients') else None
        if eft_coeffs is not None:
            # Check to see if the ordering of WCs for this sample matches what want
            if self._samples[dataset]['WCnames'] != self._wc_names_lst:
                eft_coeffs = efth.remap_coeffs(self._samples[dataset]['WCnames'], self._wc_names_lst, eft_coeffs)
        eft_w2_coeffs = efth.calc_w2_coeffs(eft_coeffs,self._dtype) if (self._do_errors and eft_coeffs is not None) else None

        # Initialize the out object
        hout = self.accumulator


        ######### Lumi Mask for Data #########        
        golden_json_path = {
            "2016": topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
            "2016APV": topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
            "2017": topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
            "2018": topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
        }
        lumi_mask = LumiMask(golden_json_path[year])(events.run,events.luminosityBlock)

        # if year == "2016" or year == "2016APV":
        #     golden_json_path = topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt")
        # elif year == "2017":
        #     golden_json_path = topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt")
        # elif year == "2018":
        #     golden_json_path = topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt")
        # else:
        #     raise ValueError(f"Error: Unknown year \"{year}\".")
        # lumi_mask = LumiMask(golden_json_path)(events.run,events.luminosityBlock)


        ######### Initialize Objects #########
        met  = events.MET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet 

        # An array of lenght events that is just 1 for each event
        events.nom = ak.ones_like(events.MET.pt)




        ######### Lepton Selection ##########
        leptonSelection = tt_os.Run2LeptonSelection()

        ele['isSelE']=leptonSelection.is_sel_ele(ele)
        mu['isSelM']=leptonSelection.is_sel_muon(mu)

        ele_good = ele[ele.isSelE]
        mu_good = mu[mu.isSelM]

        leps = ak.concatenate([ele_good, mu_good], axis=1)
        leps_sorted = leps[ak.argsort(leps.pt, axis=-1,ascending=False)] 
        leps_sorted = ak.pad_none(leps_sorted, 2)
        l0 = leps_sorted[:,0]
        l1 = leps_sorted[:,1]

        # nleps = ak.num(leps)


        ######### Systematics #########
        # wgt_correction_syst_lst = [
            # 'FSRUp', 'FSRDown', 'ISRUp', 'ISRDown', 'renormUp', 'renormDown', 'factUp', 'factDown', # Theory systs
            # ]

        # data_syst_lst = []

        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData
        # Note: add() will generally modify up/down weights, so if these are needed for any reason after this point, we should instead pass copies to add()
        # Note: Here we will add to the weights object the SFs that do not depend on any of the forthcoming loops

        # weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)

        if not isData: 
            # If this is no an eft sample, get the genWeight
            if eft_coeffs is None: 
                genw = events['genWeight']
            else: 
                genw = np.ones_like(events['event'])


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


        ######### Add variables into event object so that they persist #########
        events['njets'] = njets 
        events['leps_pt_sorted'] = leps_sorted

        # events['nbtagsm'] = nbtagsm 

        tt_es.addLepCatMasks(events)

        ######### Create Lepton Categories ##########
        select_cat_dict = None
        with open(ttbarEFT_path("params/channels.json"), "r") as ch_json_test:
            select_cat_dict = json.load(ch_json_test)

        CR_cat_dict = select_cat_dict['CR_CHANNELS_JETS'] 


        ######### Selection Masks #########
        # create trigger mask
        pass_trg = tc_es.trg_pass_no_overlap(events, isData, dataset, str(year), tt_es.triggers_dict, tt_es.exclude_triggers_dict)
        # at_least_two_leps = ak.fill_none(nleps>=2, False)

        # b jet masks
        # bmask_exactly0med = (nbtagsm==0) 
        # bmask_exactly1med = (nbtagsm==1) 
        # bmask_exactly2med = (nbtagsm==2) 
        # bmask_atleast2med = (nbtagsm>=2) 

        # Charge masks
        chargel0_p = ak.fill_none(((l0.charge)>0),False)
        chargel0_m = ak.fill_none(((l0.charge)<0),False)
        charge2l_os = ak.fill_none(((l0.charge+l1.charge)==0),False)
        charge2l_ss = ak.fill_none(((l0.charge+l1.charge)!=0),False)



        ######### Store boolean masks with PackedSelection ##########
        selections = PackedSelection(dtype='uint64')
        
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


        ######### Variables for the dense axes of the hists ##########

        ptll = (l0+l1).pt
        mll = (l0+l1).mass

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

        ######## Normalizations ########
        if not isData: 
            lumi = 1000.0*get_lumi(year)
            norm = (xsec/sow)*lumi

            if eft_coeffs is None:
                genw = events["genWeight"]
            else:
                genw = np.ones_like(events['event'])

            weights = norm*genw

        else: 
            weights = np.ones_like(events['event'])


        ########## Fill the histograms ##########

        # selections used by all categories
        # base_cuts_list = ['pass_trg', 'jetvetomap', 'bmask_exactly0med']

        # loop through categories adding selections for each
        for lep_cat in CR_cat_dict.keys():

            for jet_cat in CR_cat_dict[lep_cat]['jet_list']: 
                cuts_list = ['jetvetomap', '2los', 'bmask_exactly0med']

                if isData:
                    cuts_list.append('is_good_lumi')
                
                cuts_list.append(lep_cat)
                cuts_list.append(jet_cat)

                ch_name = lep_cat

                event_selection_mask = selections.all(*(cuts_list))
                weights_cut = weights[event_selection_mask]
                eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

                for dense_axis_name, dense_axis_vals in dense_axis_variables.items():
                    # if the category requires zero bjets, don't fill jet histograms
                    if (jet_cat == 'exactly_0j') and (dense_axis_name in jet_variables):
                        continue
                    if dense_axis_name not in self._hist_lst:
                        print(f"Skipping \"{dense_axis_name}\", it is not in the list of hists to include")
                        continue 

                    # Fill the histos
                    axes_fill_info_dict = {
                        dense_axis_name : dense_axis_vals[event_selection_mask],
                        "channel"       : ch_name,
                        "process"       : histAxisName,
                        "weight"        : weights_cut,
                        "eft_coeff"     : eft_coeffs_cut,
                    }

                    hout[dense_axis_name].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator


