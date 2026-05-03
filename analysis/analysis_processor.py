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
from ttbarEFT.modules.corrections import AttachElectronSF, AttachMuonSF, AttachMuonTrigSF, ApplyMuonPtCorr, ApplyJetVetoMaps

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
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)
        # appl_axis = hist.axis.StrCategory([], name="appl", label=r"AR/SR", growth=True)

        # Set the booleans
        self._do_errors = do_errors # Whether to calculate and store the w**2 coefficients
        self._do_systematics = do_systs # Whether to process systematic samples
        self._split_by_lepton_flavor = split_by_lepton_flavor # Whether to keep track of lepton flavors individually
        self._skip_signal_regions = skip_signal_regions # Whether to skip the SR categories
        self._skip_control_regions = skip_control_regions # Whether to skip the CR categories

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
                chan_axis,
                syst_axis,
                dense_axis,
                wc_names = wc_names_lst, 
                label=r'Events',
            )

            histograms[name+'_sumw2'] = HistEFT(
                proc_axis, 
                chan_axis,
                syst_axis,
                sumw2_axis,
                wc_names = wc_names_lst, 
                label=r'Events',
            )

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
        histAxisName    = self._samples[dataset]['histAxisName']
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']

        isEFT = hasattr(events, 'EFTfitCoefficients')    
        assert not (isEFT and isData), f"isEFT and isData cannot both be True. Check input samples."


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


        ######### Initialize Objects #########
        met  = events.MET
        ele  = events.Electron
        mu   = events.Muon
        tau  = events.Tau
        jets = events.Jet 

        leptonSelection = tt_os.Run2LeptonSelection()

        # An array of length events that is just 1 for each event
        events.nom = ak.ones_like(events.MET.pt)


        ######### Electron Selection ##########
        
        ele['isSelE']=leptonSelection.is_sel_ele(ele)
        ele_good = ele[ele.isSelE]

        AttachElectronSF(ele_good, year) 


        ######### Muon Selection ##########
        
        mu['pt'] = ApplyMuonPtCorr(mu, year, isData)
        mu['isSelM']=leptonSelection.is_sel_muon(mu)
        
        mu_good = mu[mu.isSelM]
        AttachMuonSF(mu_good, year)


        ######### Lepton Selection ##########
        leps = ak.concatenate([ele_good, mu_good], axis=1)
        leps_sorted = leps[ak.argsort(leps.pt, axis=-1,ascending=False)] 

        events['leps_pt_sorted'] = leps_sorted


        ######### Systematics #########

        obj_correction_syst_lst = [] #TO-DO

        wgt_correction_syst_lst = [
            'SF_ele_up', 'SF_ele_down', 'SF_muon_up', 'SF_muon_down' # Exp systs
            # 'FSRUp', 'FSRDown', 'ISRUp', 'ISRDown', 'renormUp', 'renormDown', 'factUp', 'factDown', # Theory systs #TO-DO
            ]

        data_syst_lst = []

        tt_es.addLepSFs(events)

        # These weights can go outside of the outside sys loop since they do not depend on pt of mu or jets
        # We only calculate these values if not isData
        # Note: add() will generally modify up/down weights, so if these are needed for any reason after this point, we should instead pass copies to add()
        # Note: Here we will add to the weights object the SFs that do not depend on any of the forthcoming loops

        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)

        if not isData: 
            # If this is no an eft sample, get the genWeight
            if eft_coeffs is None: 
                genw = events['genWeight']
            else: 
                genw = np.ones_like(events['event'])

            lumi = 1000.0*get_lumi(year)
            norm = genw*(xsec/sow)*lumi
            weights_obj_base.add('norm', norm)

            # AttachPSWeights(events)
            # AttachScaleWeights(events) 
            # AttachPdfWeights(events)

            # weights_obj_base.add('ISR', events.nom, events.ISRUp*(sow/sow_ISRUp), events.ISRDown*(sow/sow_ISRDown))
            # weights_obj_base.add('FSR', events.nom, events.FSRUp*(sow/sow_FSRUp), events.FSRDown*(sow/sow_FSRDown))
            # # renorm/fact scale  -- corrections come from AttachScaleWeights
            # weights_obj_base.add('renorm', events.nom, events.renormUp*(sow/sow_renormUp), events.renormDown*(sow/sow_renormDown))
            # weights_obj_base.add('fact', events.nom, events.factUp*(sow/sow_factUp), events.factDown*(sow/sow_factDown))
            # # Prefiring and PU (note prefire weights only available in nanoAODv9 and for Run2)
            # weights_obj_base.add('PreFiring', *l1prefiring_args) #Run3 ready
            # weights_obj_base.add('PU', tc_cor.GetPUSF((events.Pileup.nTrueInt), year), tc_cor.GetPUSF(events.Pileup.nTrueInt, year, 'up'), tc_cor.GetPUSF(events.Pileup.nTrueInt, year, 'down')) #Run3 ready

            # current understanding is that lepSFs are applied to all event channels, and don't affect any kinematics so they can go here safely
            weights_obj_base.add('lepSF_ele', events.SF_2l_ele, copy.deepcopy(events.SF_2l_ele_up), copy.deepcopy(events.SF_2l_ele_down))
            weights_obj_base.add('lepSF_muon', events.SF_2l_muon, copy.deepcopy(events.SF_2l_muon_up), copy.deepcopy(events.SF_2l_muon_down)) 

            # should be able to add triggerSF here also 

        ######### Create Lepton Categories ##########
        select_cat_dict = None
        with open(ttbarEFT_path("params/channels.json"), "r") as ch_json_test:
            select_cat_dict = json.load(ch_json_test)

        CR_cat_dict = select_cat_dict['CR_CHANNELS_JETS'] 

 
        ######### The rest of the processor is inside this loop over systs that affect object kinematics  ###########

        # If we're doing systematics and this isn't data, we will loop over the obj_correction_syst_lst list
        if self._do_systematics and not isData: syst_var_list = ['nominal'] + obj_correction_syst_lst
        # Otherwise loop juse once, for nominal
        else: syst_var_list = ['nominal'] 

        for syst_var in syst_var_list:
            # Make a copy of the base weights object, so that each time through the loop we do not double count systs
            # In this loop over systs that impact kinematics, we will add to the weights objects the SFs that depend on the object kinematics
            weights_obj_base_for_kinematic_syst = copy.deepcopy(weights_obj_base)


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
            # events['nbtagsm'] = nbtagsm 

            tt_es.addLepCatMasks(events)


            ######### Create objects for dense axes ##########
            leps_sorted = ak.pad_none(leps_sorted, 2)
            l0 = leps_sorted[:,0]
            l1 = leps_sorted[:,1]

            ptll = (l0+l1).pt
            mll = (l0+l1).mass


            ######### Event weights that do not depend on the lep cat ##########
            # if not isData: 

                # btag SFs and systematics go here
               

            # weights_dict = {}
            # for ch_name in lep_cats: 
            #     weights_dict[ch_name] = copy.deepcopy(weights_obj_base_for_kinematic_syst)


            ######### Selection Masks #########
            # pass_trg = tt_es.trg_pass_no_overlap(events, isData, dataset, str(year), tt_es.triggers_dict, tt_es.exclude_triggers_dict, lep_cat)
            pass_trg = tc_es.trg_pass_no_overlap(events, isData, dataset, str(year), tt_es.triggers_dict, tt_es.exclude_triggers_dict)
            # at_least_two_leps = ak.fill_none(nleps>=2, False)

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

            # setup list of systematics to loop through 
            wgt_var_lst = ['nominal']
            if self._do_systematics: 
                if not isData: 
                    if (syst_var != 'nominal'): 
                            # In this case, we are dealing with systs that change the kinematics of the objs (e.g. JES)
                            # So we don't want to loop over up/down weight variations here
                            wgt_var_lst = [syst_var]
                    else: 
                        # Otherwise we want to loop over the up/down weight variations
                        wgt_var_lst = wgt_var_lst + wgt_correction_syst_lst 

                else: 
                    # This is data, so we want to loop over just up/down variations relevant for data (eg FF)
                    wgt_var_lst = wgt_var_lst + data_syst_lst

            # Loop through systematics and fill histograms
            for wgt_fluct in wgt_var_lst:

                # currently, weights_obj is the same for all channels
                # later replace this with weights[ch_name] inside of a loop over channels if needed 
                weights_object = weights_obj_base_for_kinematic_syst 
                if (wgt_fluct == "nominal") or (wgt_fluct in obj_correction_syst_lst):
                    # In the case of "nominal", or the jet energy systematics, no weight systematic variation is used
                    weight = weights_object.weight(None)  
                else: 
                    # Get the weight variation from the Weights object 
                    if wgt_fluct in weights_object.variations:
                        weight = weights_object.weight(wgt_fluct)
                    else: 
                        # If there is no up/down in this weight category, we don't want to fill a histogram 
                        continue 

                # This is a check ot make sure we guard against any unintentional variations being applied to data
                if self._do_systematics and isData:
                    # the up/down variations should correspond to only the ones in the data list 
                    if weights_object.variations != set(data_syst_lst): 
                        raise Exception(f"Error: Unexpected wgt variations for data! Expected \"{set(data_syst_lst)}\" but have \"{weights_object.variations}\".")

                # loop through categories adding selections for each
                for lep_cat in CR_cat_dict.keys():

                    for jet_cat in CR_cat_dict[lep_cat]['jet_list']: 
                        # masks that are applied to all categories
                        cuts_list = ['jetvetomap', '2los', 'bmask_exactly0med']

                        if isData:
                            cuts_list.append('is_good_lumi')
                        
                        cuts_list.append(lep_cat)
                        cuts_list.append(jet_cat)

                        # ch_name = f"{lep_cat}_{jet_cat}"
                        ch_name = f"{lep_cat}"

                        event_selection_mask = selections.all(*(cuts_list))
                        weights_cut = weight[event_selection_mask]
                        eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

                        for dense_axis_name, dense_axis_vals in dense_axis_variables.items():
                            # if the category requires zero jets, don't fill jet histograms
                            if (jet_cat == 'exactly_0j') and (dense_axis_name in jet_variables):
                                print(f"Skipping '{dense_axis_name}' in category '{ch_name}_{jet_cat}'. Jet histograms are not filled for categories that don't require a jet")
                                continue

                            if dense_axis_name not in self._hist_lst:
                                print(f"Skipping \"{dense_axis_name}\", it is not in the list of hists to include")
                                continue                         

                            # Fill the histos
                            axes_fill_info_dict = {
                                dense_axis_name : dense_axis_vals[event_selection_mask],
                                "channel"       : ch_name,
                                "process"       : histAxisName,
                                "systematic"    : wgt_fluct,
                                "weight"        : weights_cut,
                                "eft_coeff"     : eft_coeffs_cut,
                            }

                            hout[dense_axis_name].fill(**axes_fill_info_dict)

                            if self._do_errors: 
                                axes_fill_info_dict = {
                                    dense_axis_name+"_sumw2" : dense_axis_vals[event_selection_mask],
                                    "channel"       : ch_name,
                                    "process"       : histAxisName,
                                    "systematic"    : wgt_fluct,
                                    "weight"        : np.square(weights_cut),
                                    "eft_coeff"     : None,
                                }

                                hout[dense_axis_name+"_sumw2"].fill(**axes_fill_info_dict)

        return hout

    def postprocess(self, accumulator):
        return accumulator


