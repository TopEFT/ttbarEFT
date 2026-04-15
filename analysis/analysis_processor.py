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
import ttbarEFT.modules.corrections as tt_cor 
from ttbarEFT.modules.processor_tools import calc_eft_weights
from ttbarEFT.modules.processor_tools import get_syst_lists

from topcoffea.modules.get_param_from_jsons import GetParam

get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_tt_param = GetParam(ttbarEFT_path("params/params.json"))

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')

# compare against topcoffea.modules.eft_helper.calc_eft_weights
# def calc_eft_weights(eft_coeffs, wc_vals):
#     '''
#     Returns an array that contains the event weight for each event.
#     eft_coeffs: Array of eft fit coefficients for each event
#     wc_vals: wilson coefficient values desired for the event weight calculation, listed in the same order as the wc_lst
#              such that the multiplication with eft_coeffs is correct
#              The correct ordering can be achieved with the order_wc_values function
#     '''
#     event_weight = np.empty_like(eft_coeffs)

#     wcs = np.hstack((np.ones(1),wc_vals))
#     wc_cross_terms = []
#     index = 0
#     for j in range(len(wcs)):
#         for k in range (j+1):
#             term = wcs[j]*wcs[k]
#             wc_cross_terms.append(term)
#     event_weight = np.sum(np.multiply(wc_cross_terms, eft_coeffs), axis=1)

#     return event_weight


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, lep_cat, wc_names_lst=[], hist_lst=None, do_errors=False, syst_list=[], dtype=np.float32):
        self._samples = samples
        self._lep_cat = lep_cat
        self._wc_names_lst = wc_names_lst
        self._do_errors = do_errors

        self._syst_list = syst_list
        # self._do_systematics = syst_list is not None # only do systematics if sys_list is not None
        self._dtype = dtype 
        

        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)

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

        if not isData: 
            # if 'sow_ISRUp' in self._samples[dataset].keys(): 
            sow_ISRUp          = self._samples[dataset]["nSumOfWeights_ISRUp"          ]
            sow_ISRDown        = self._samples[dataset]["nSumOfWeights_ISRDown"        ]
            sow_FSRUp          = self._samples[dataset]["nSumOfWeights_FSRUp"          ]
            sow_FSRDown        = self._samples[dataset]["nSumOfWeights_FSRDown"        ]
            sow_renormUp       = self._samples[dataset]["nSumOfWeights_renormUp"       ]
            sow_renormDown     = self._samples[dataset]["nSumOfWeights_renormDown"     ]
            sow_factUp         = self._samples[dataset]["nSumOfWeights_factUp"         ]
            sow_factDown       = self._samples[dataset]["nSumOfWeights_factDown"       ]
            sow_renormfactUp   = self._samples[dataset]["nSumOfWeights_renormfactUp"   ]
            sow_renormfactDown = self._samples[dataset]["nSumOfWeights_renormfactDown" ]

        print(f"\n\n")
        print(f"histAxisName: {histAxisName}")
        print(f"dataset: {dataset}")
        print(f"year: {year}")
        print(f"xsec: {xsec}")
        print(f"\n\n")

        isEFT = hasattr(events, 'EFTfitCoefficients')    
        assert not (isEFT and isData), f"isEFT and isData cannot both be True. Check input samples."

        lep_cat = self._lep_cat

        run_era = None
        datasets = ["Muon", "SingleMuon", "SingleElectron", "EGamma", "MuonEG", "DoubleMuon", "DoubleElectron", "DoubleEG"]
        for d in datasets:
            if dataset.startswith(d):
                dataset_string = dataset.split('_')
                dataset = dataset_string[0]
                run_era = dataset_string[1]


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

        ######### Data Selections #########
        if isData:
            # Lumi Mask for Data    
            golden_json_path = {
                "2016"      : topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
                "2016APV"   : topcoffea_path("data/goldenJsons/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"),
                "2017"      : topcoffea_path("data/goldenJsons/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"),
                "2018"      : topcoffea_path("data/goldenJsons/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"),
            }
            lumi_mask = LumiMask(golden_json_path[year])(events.run,events.luminosityBlock) 


        obj_correction_syst_lst = tt_cor.get_supported_jet_systematics(year, isData=isData, era=run_era)
        event_weight_variations, kinematic_variations = get_syst_lists(year=year, isData=isData, syst_list=self._syst_list, run_era=None)

        print(f"\n\n")
        print(f"list of systematics to run over: \n\tevent_weight_variations = {event_weight_variations}, \n\tkinematic_variations = {kinematic_variations}")
        print(f"\n\n")

        ######### Load Event Categories ##########
        cat_dict = None
        with open(ttbarEFT_path("params/channels.yaml"), "r") as f:
            cat_dict=yaml.safe_load(f)

        # CR_cat_dict = cat_dict['CR_CHANNELS']
        CR_cat_dict = cat_dict['CR_CHANNELS_allb']
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

        ######### Muon Selection ##########
        mu['pt'] = tt_cor.ApplyMuonPtCorr(mu, year, isData)
        mu['isGoodMuon']=leptonSelection.is_sel_muon(mu)
        mu_good = mu[mu.isGoodMuon]

        ######### Lepton Selection ##########
        # add lepton scale factors and trigger efficiencies
        if not isData: 
            tt_cor.AttachElectronSF(ele_good, year)    
            tt_cor.AttachMuonSF(mu_good, year)    
            tt_cor.AttachElecTrigEff(ele_good, year)
            tt_cor.AttachMuonTrigEff(mu_good, year) 

        leps = ak.concatenate([ele_good, mu_good], axis=1)
        leps_sorted = leps[ak.argsort(leps.pt, axis=-1,ascending=False)] 


        ######### Jet Selections #########
        jets['isClean'] = tt_os.isClean(jets, ele_good, drmin=0.4)& tt_os.isClean(jets, mu_good, drmin=0.4)
        cleanedJets = jets[jets.isClean]
        cleanedJets = ak.fill_none(cleanedJets, [])

        # Medium DeepJet WP
        medium_tag = "btag_wp_medium_" + year.replace("201", "UL1")
        btagwpm = get_tc_param(medium_tag)


        ######### Add variables to EVENTS #########
        events['leps_pt_sorted'] = leps_sorted
        tt_es.addLepCatMasks(events) 
        tt_es.add2losMask(events, year, isData)

        ######## Create objects for dense axes ##########
        leps_sorted = ak.pad_none(leps_sorted, 2)
        l0 = leps_sorted[:,0]
        l1 = leps_sorted[:,1]

        ptll = (l0+l1).pt
        mll = (l0+l1).mass


        ######### Selection Masks that aren't dependent on object corrections #########
        pass_trg = tt_es.trg_pass_no_overlap(events, isData, dataset, str(year), tt_es.triggers_dict, tt_es.exclude_triggers_dict, lep_cat)

        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)

        if not isData: 
            if eft_coeffs is None:                      
                genw = events['genWeight']
            else:                                       # If this is not an eft sample, get the genWeight
                genw = np.ones_like(events['event'])

            lumi = 1000.0*get_lumi(year)
            norm = genw*(xsec/sow)*lumi
            weights_obj_base.add('norm', norm)

            weights_obj_base.add('elecID', *tt_cor.Get_ElecIDSF(events))
            weights_obj_base.add('muonID', *tt_cor.Get_MuonIDSF(events))
            weights_obj_base.add('muonISO', *tt_cor.Get_MuonISOSF(events))

            weights_obj_base.add('trigSF', *tt_cor.GetTrigSF(events, lep_cat)) # a bit misleading, ee is trigger efficiencies so up/down is set to ones

            weights_obj_base.add('L1prefire', events.L1PreFiringWeight.Nom, events.L1PreFiringWeight.Up, events.L1PreFiringWeight.Dn)
            weights_obj_base.add('PU', tt_cor.GetPUSF(events.Pileup.nTrueInt, year), tt_cor.GetPUSF(events.Pileup.nTrueInt, year, 'up'), tt_cor.GetPUSF(events.Pileup.nTrueInt, year, 'down'))

            tt_cor.AttachScaleWeights(events) # Matrix Element uncertainties
            weights_obj_base.add('renorm', events.nom, events.renormUp*(sow/sow_renormUp), events.renormDown*(sow/sow_renormDown))
            weights_obj_base.add('fact', events.nom, events.factUp*(sow/sow_factUp), events.factDown*(sow/sow_factDown))
            
            tc_cor.AttachPSWeights(events) #QScale uncertainties
            weights_obj_base.add('ISR', events.nom, events.ISRUp*(sow/sow_ISRUp), events.ISRDown*(sow/sow_ISRDown))
            weights_obj_base.add('FSR', events.nom, events.FSRUp*(sow/sow_FSRUp), events.FSRDown*(sow/sow_FSRDown))


        # for Run2, Jet Corrections are applied to Data, only run this on MC
        if not isData:
            # Medium DeepJet WP lookup functions
            btag_eff_lookup_m = tt_cor.GetBtagEffLookup(year, wp='medium')
            light_btag_SF_lookup = tt_cor.GetBtagSFLookup(wp='M',year=year, method='deepJet_incl')
            bc_btag_SF_lookup = tt_cor.GetBtagSFLookup(wp='M',year=year, method='deepJet_comb')

            raw_met = met
            cleanedJets['pt_orig'] = cleanedJets.pt     # NECESSARY FOR MET CORRECTIONS LATER 
            # cleanedJets["pt_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.pt
            # cleanedJets["mass_raw"] = (1 - cleanedJets.rawFactor)*cleanedJets.mass
            # # cleanedJets["rho"] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, cleanedJets.pt)[0] #THIS LINE BREAKS THE JETS BUT NOT IN A WAY THAT FAILS
            rho_jagged = ak.ones_like(cleanedJets.pt) * events.fixedGridRhoFastjetAll
            cleanedJets = ak.with_field(cleanedJets, rho_jagged, "Rho")
            cleanedJets["pt_gen"] = ak.values_astype(ak.fill_none(cleanedJets.matched_gen.pt, 0), np.float32)
            cleanedJets = tt_cor.ApplyJetCorrections(year, corr_type='jets', isData=isData, era=run_era).build(cleanedJets)
            corrected_met = tt_cor.ApplyJetCorrections(year, corr_type='met', isData=isData, era=run_era).build(MET=raw_met, corrected_jets=cleanedJets)

        ######### The rest of the processor is inside this loop over systs that affect object kinematics  ###########
        print(f"\n\n")
        print(f"kinematic_variations: {kinematic_variations}")
        print(f"")
        print(f"\n\n")

        for kinematic_var in kinematic_variations: 

            if isData: 
                cleanedJets['isGood'] = tt_os.is_pres_jet(cleanedJets)
                goodJets =  cleanedJets[cleanedJets.isGood]

                jet_veto_map = tt_cor.ApplyJetVetoMaps(goodJets, year)

                njets = ak.num(goodJets)
                jets_sorted = goodJets[ak.argsort(goodJets.pt, axis=-1,ascending=False)]
                jets_sorted = ak.pad_none(jets_sorted, 1)
                j0 = jets_sorted[:,0]
                # ht = ak.sum(goodJets.pt,axis=-1)
                isBtagJetsMedium = (goodJets.btagDeepFlavB > btagwpm)
                nbtagsm = ak.num(goodJets[isBtagJetsMedium])


            elif not isData:

                weights_obj_base_for_kinematic_syst = copy.deepcopy(weights_obj_base)

                if kinematic_var == 'nominal': 
                    cleanedJets['isGood'] = tt_os.is_pres_jet(cleanedJets)
                    goodJets =  cleanedJets[cleanedJets.isGood]
                    met = corrected_met

                else: 
                    correctedJets = tt_cor.ApplyJetSystematics(year=year, corr_type='jets', original_obj=cleanedJets, syst_var=kinematic_var)
                    met = tt_cor.ApplyJetSystematics(year=year, corr_type='met', original_obj=corrected_met, syst_var=kinematic_var)

                    correctedJets['isGood'] = tt_os.is_pres_jet(correctedJets)
                    goodJets = correctedJets[correctedJets.isGood]

                weights_obj_base_for_kinematic_syst.add('jetPuID', tt_cor.GetJetPuIDSF(year, goodJets, var='nom'), tt_cor.GetJetPuIDSF(year, goodJets, var='up'), tt_cor.GetJetPuIDSF(year, goodJets, var='down'))
                jet_veto_map = tt_cor.ApplyJetVetoMaps(goodJets, year)    

                njets = ak.num(goodJets)
                jets_sorted = goodJets[ak.argsort(goodJets.pt, axis=-1,ascending=False)]
                jets_sorted = ak.pad_none(jets_sorted, 1)
                j0 = jets_sorted[:,0]                                       # ht = ak.sum(goodJets.pt,axis=-1)
                isBtagJetsMedium = (goodJets.btagDeepFlavB > btagwpm)
                nbtagsm = ak.num(goodJets[isBtagJetsMedium])

                # nominal btag SF 
                light_jets = goodJets[goodJets.hadronFlavour == 0]
                light_eff = btag_eff_lookup_m(light_jets)
                light_bSF = light_btag_SF_lookup(jet_collection=light_jets, syst='central')
                light_btagweight = tt_cor.GetBtag_method1a_wgt_singlewp(light_eff, light_bSF, passes_tag=light_jets.btagDeepFlavB > btagwpm)

                bc_jets = goodJets[goodJets.hadronFlavour > 0]
                bc_eff = btag_eff_lookup_m(bc_jets)
                bc_bSF = bc_btag_SF_lookup(jet_collection=bc_jets, syst='central')
                bc_btagweight = tt_cor.GetBtag_method1a_wgt_singlewp(bc_eff, bc_bSF, passes_tag=bc_jets.btagDeepFlavB > btagwpm)

                btag_eventweight = light_btagweight*bc_btagweight
                weights_obj_base_for_kinematic_syst.add('btagSF', btag_eventweight)


                if (kinematic_var=='nominal') and ('btagSFbc_correlatedUp' in event_weight_variations): 
                    print(f"running over btagSF variations")

                    for b_syst in ['bc_correlated', 'light_correlated', f"bc_{year}", f"light_{year}"]:
                        if b_syst.endswith("correlated"): 
                            corrtype = "correlated"
                        else: 
                            corrtype = "uncorrelated"

                        if b_syst.startswith("light_"):
                            sf_lookup = light_btag_SF_lookup
                            syst_jets = light_jets
                            btag_effM = light_eff
                            fixed_btag_w = bc_btagweight
                        elif b_syst.startswith("bc_"):
                            sf_lookup = bc_btag_SF_lookup
                            syst_jets = bc_jets
                            btag_effM = bc_eff
                            fixed_btag_w = light_btagweight


                        btag_SF_up = sf_lookup(jet_collection=syst_jets,syst=f'up_{corrtype}')
                        btag_SF_down = sf_lookup(jet_collection=syst_jets,syst=f'down_{corrtype}')

                        btagweight_up = tt_cor.GetBtag_method1a_wgt_singlewp(btag_effM, btag_SF_up, passes_tag=syst_jets.btagDeepFlavB > btagwpm)
                        btagweight_down = tt_cor.GetBtag_method1a_wgt_singlewp(btag_effM, btag_SF_down, passes_tag=syst_jets.btagDeepFlavB > btagwpm)

                        event_weight_up = fixed_btag_w * btagweight_up
                        event_weight_down = fixed_btag_w * btagweight_down

                        # TODO: maybe change events.nom to btag_eventweight? depends on how this will be used in combine
                        weights_obj_base_for_kinematic_syst.add(f'btagSF{b_syst}', events.nom, event_weight_up/btag_eventweight, event_weight_down/btag_eventweight)

            # add HEM veto using "good" leptons and jets
            HEM_veto_mask, HEM_event_weight = tt_es.getHemMask(events, year, isData, jets=goodJets)
            
            if not isData:
                weights_obj_base_for_kinematic_syst.add(f'HEM', HEM_event_weight)
            
            ######### Store boolean masks with PackedSelection ##########
            selections = PackedSelection(dtype='uint64')

            if isData: 
                selections.add('is_good_lumi', lumi_mask)

            selections.add('pass_trg', pass_trg)
            # selections.add('2los', events.is2los)

            if isData: 
                selections.add("ee",  (events.is_ee & events.is2los & pass_trg))    # data always has trig pass requirement
            else: 
                selections.add("ee",  (events.is_ee & events.is2los))               # MC for ee has efficiencies, so no pass_trg requirement
                
            selections.add("em",  (events.is_em & events.is2los & pass_trg))        # MC for emu has SF, so use pass_trg requirement 
            selections.add("mm",  (events.is_mm & events.is2los & pass_trg))        # MC for mumu has SF, so use pass_trg requirement

            selections.add('jetvetomap', (jet_veto_map == 0))
            selections.add('HEMvetomap', HEM_veto_mask)

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
            dense_axis_variables['mllZ'] = mll
            dense_axis_variables['ptll'] = ptll
            dense_axis_variables['l0pt'] = l0.pt
            dense_axis_variables['l0eta'] = l0.eta 
            dense_axis_variables['l0phi'] = l0.phi 
            dense_axis_variables['l1pt'] = l1.pt
            dense_axis_variables['l1eta'] = l1.eta 
            dense_axis_variables['l1phi'] = l1.phi
            dense_axis_variables['j0pt'] = j0.pt
            dense_axis_variables['j0eta'] = j0.eta
            dense_axis_variables['j0phi'] = j0.phi
            dense_axis_variables['MET'] = met.pt

            jet_variables = ['j0pt', 'j0eta', 'j0phi']


            ########## Fill the histograms ##########
            wgt_var_lst = ["nominal"]
            if self._syst_list and not isData: #if systematics were provided and it's MC create list of variations 
                if (kinematic_var != "nominal"):
                    # in this case, we are dealing with systs that change the kinematics of objects
                    # we don't want to loop over up/down weight variations here
                    wgt_var_lst = [kinematic_var]
                else: 
                    # in this case we want to loop over the up/down event weight variations
                    wgt_var_lst = wgt_var_lst + event_weight_variations


            lep_cat_channels = CR_cat_dict[lep_cat]
            for wgt_fluct in wgt_var_lst: 

                if isData:
                    weight = np.ones_like(events.event)
                elif (wgt_fluct == "nominal") or (wgt_fluct in obj_correction_syst_lst):
                    weight = weights_obj_base_for_kinematic_syst.weight(None) 
                else: 
                    if wgt_fluct in weights_obj_base_for_kinematic_syst.variations:
                        weight = weights_obj_base_for_kinematic_syst.weight(wgt_fluct)
                    else: 
                        continue        # if there is no up/down fluctuation for this category, don't fill a hist

                for chan_id, chan_settings in lep_cat_channels.items():
                    chan_name = chan_settings['name']
                    mask_list = chan_settings['masks']

                # for jet_cat in CR_cat_dict[lep_cat]['jet_list']: 
                    # masks that are applied to all categories
                    cuts_list = ['jetvetomap', 'HEMvetomap']

                    if isData:
                        cuts_list.append('is_good_lumi')
                    
                    cuts_list.append(lep_cat)
                    cuts_list.extend(mask_list)

                    print(f"Filling bin '{chan_name}' using masks: {cuts_list}")

                    event_selection_mask = selections.all(*(cuts_list))
                    eft_coeffs_cut = eft_coeffs[event_selection_mask] if eft_coeffs is not None else None

                    ### PDF Weights
                    # pdfWeights = events.LHEPdfWeight
                    #     for i in range(100): #maybe 101 instead of 100? 
                    #     # weight = orig_weight[event_selection_mask] * pdfWeights[event_selection_mask]

                    # channel_name = lep_cat

                    for dense_axis_name, dense_axis_vals in dense_axis_variables.items():
                        # if the category requires zero jets, don't fill jet histograms
                        if ('exactly_0j' in cuts_list) and (dense_axis_name in jet_variables):
                            # print(f"Skipping '{dense_axis_name}' in category '{lep_cat}_{jet_cat}'. Jet histograms are not filled for categories that don't require a jet")
                            continue

                        if dense_axis_name not in self._hist_lst:
                            print(f"Skipping \"{dense_axis_name}\", it is not in the list of hists to include")
                            continue                       

                        # Fill the histos
                        axes_fill_info_dict = {
                            dense_axis_name : dense_axis_vals[event_selection_mask],
                            "process"       : histAxisName,
                            'channel'       : chan_name,
                            "systematic"    : wgt_fluct,
                            "weight"        : weight[event_selection_mask],
                            "eft_coeff"     : eft_coeffs_cut,
                        }

                        hout[dense_axis_name].fill(**axes_fill_info_dict)

                        if self._do_errors and wgt_fluct=='nominal': 

                            if eft_coeffs is not None:
                                event_weights_SM = calc_eft_weights(eft_coeffs,np.zeros(len(self._wc_names_lst)))
                                sumw2 = np.square(weight*event_weights_SM)
                            else: 
                                sumw2 = np.square(weight)

                            sumw2axes_fill_info_dict = {
                                dense_axis_name             : dense_axis_vals[event_selection_mask],
                                'process'                   : histAxisName,
                                'channel'                   : chan_name,
                                'systematic'                : 'sumw2',
                                'weight'                    : sumw2[event_selection_mask],
                                'eft_coeff'                 : None,
                            }

                            hout[dense_axis_name].fill(**sumw2axes_fill_info_dict)

        return hout
        

    def postprocess(self, accumulator):
        return accumulator

