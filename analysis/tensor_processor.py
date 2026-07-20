#!/usr/bin/env python
import copy
import numpy as np
import awkward as ak

import coffea, os, random, torch, time

# from mt2 import mt2

from coffea import processor
from coffea.analysis_tools import PackedSelection
from coffea.lumi_tools import LumiMask

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoAODSchema

from topcoffea.modules.paths import topcoffea_path
import topcoffea.modules.eft_helper as efth
import topcoffea.modules.corrections as tc_cor
import topcoffea.modules.event_selection as tc_es

from ttbarEFT.modules.paths import ttbarEFT_path
from ttbarEFT.modules.analysis_tools import get_lumi
import ttbarEFT.modules.object_selection as tt_os
import ttbarEFT.modules.event_selection as tt_es
import ttbarEFT.modules.corrections as tt_cor 
from ttbarEFT.modules.processor_tools import get_syst_lists
from ttbarEFT.modules.analysis_tools import TensorAccumulator

from topcoffea.modules.get_param_from_jsons import GetParam

get_tc_param = GetParam(topcoffea_path("params/params.json"))
get_tt_param = GetParam(ttbarEFT_path("params/params.json"))

NanoAODSchema.warn_missing_crossrefs = False
np.seterr(divide='ignore', invalid='ignore', over='ignore')


class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, lep_cat, outname, wc_names_lst=[], do_errors=False, doPDF=False, doSR=False, syst_list=[], dtype=np.float32):
        self._samples = samples
        self._lep_cat = lep_cat
        self._wc_names_lst = wc_names_lst
        self._do_errors = do_errors
        self._doPDF = doPDF

        self._outname = outname

        self._doCR = not doSR
        self._doSR = doSR

        self._syst_list = syst_list
        self._dtype = dtype 

        os.makedirs(f'{outname}/training',   mode=0o755, exist_ok=True)
        os.makedirs(f'{outname}/validation', mode=0o755, exist_ok=True)

    @property
    def accumulator(self):
        return self._accumulator

    @property
    def columns(self):
        return self._columns

    def process(self, events):

        out = TensorAccumulator(torch.tensor([]))

        # Dataset parameters
        dataset         = events.metadata['dataset']
        isData          = self._samples[dataset]['isData']
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']

        if not isData: 
            # if 'sow_ISRUp' in self._samples[dataset].keys(): 
            sow_ISRUp           = self._samples[dataset]["nSumOfWeights_ISRUp"          ]
            sow_ISRDown         = self._samples[dataset]["nSumOfWeights_ISRDown"        ]
            sow_FSRUp           = self._samples[dataset]["nSumOfWeights_FSRUp"          ]
            sow_FSRDown         = self._samples[dataset]["nSumOfWeights_FSRDown"        ]
            sow_renormUp        = self._samples[dataset]["nSumOfWeights_renormUp"       ]
            sow_renormDown      = self._samples[dataset]["nSumOfWeights_renormDown"     ]
            sow_factUp          = self._samples[dataset]["nSumOfWeights_factUp"         ]
            sow_factDown        = self._samples[dataset]["nSumOfWeights_factDown"       ]
            sow_renormfactUp    = self._samples[dataset]["nSumOfWeights_renormfactUp"   ]
            sow_renormfactDown  = self._samples[dataset]["nSumOfWeights_renormfactDown" ]
            sow_hdampUp         = self._samples[dataset]["nSumOfWeights_hdampUp"        ]
            sow_hdampDown       = self._samples[dataset]["nSumOfWeights_hdampDown"      ]
            sow_toppt           = self._samples[dataset]["nSumOfWeights_toppt"          ]

        print(f"\n\n")
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
        obj_correction_syst_lst += ['METunclustUp', 'METunclustDown']
        event_weight_variations, kinematic_variations = get_syst_lists(year=year, isData=isData, syst_list=self._syst_list, run_era=None)

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
            ele_good = tt_cor.AttachElectronSF(ele_good, year)    
            mu_good = tt_cor.AttachMuonSF(mu_good, year)    
            ele_good = tt_cor.AttachElecTrigEff(ele_good, year)
            mu_good = tt_cor.AttachMuonTrigEff(mu_good, year) 

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
        tt_es.addmllMasks(events)

        ######## Create objects for dense axes ##########
        leps_sorted = ak.pad_none(leps_sorted, 2)
        pMask = leps_sorted.pdgId < 0
        nMask = leps_sorted.pdgId > 0
        l0 = leps_sorted[:,0]
        l1 = leps_sorted[:,1]
        lp = leps_sorted[pMask, 0]
        ln = leps_sorted[nMask, 0]


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

            # weights_obj_base.add('elecID', *tt_cor.Get_ElecIDSF(events))
            weights_obj_base.add('elecHEEP', *tt_cor.Get_ElecHEEPSF(events))
            weights_obj_base.add('elecRECO', *tt_cor.Get_ElecRECOSF(events))
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

            weights_obj_base.add('hdamp', events.nom, (tt_cor.GetHdampReweight(events, dataset, var='up')*(sow/sow_hdampUp)), (tt_cor.GetHdampReweight(events, dataset, var='down')*(sow/sow_hdampDown)))
            
            LOtoNLO_weights = tt_cor.GetNLO_Weight(events, dataset)
            NLOtoNNLO_weights = tt_cor.GetNNLO_EventWeight(events, dataset)
            weights_obj_base.add('ttbar_toppt', LOtoNLO_weights*NLOtoNNLO_weights*(sow/sow_toppt))

        # for Run2, Jet Corrections are applied to Data, only run this on MC
        if not isData:
            # Medium DeepJet WP lookup functions
            btag_eff_lookup_m = tt_cor.GetBtagEffLookup(year, wp='medium')
            light_btag_SF_lookup = tt_cor.GetBtagSFLookup(wp='M',year=year, method='deepJet_incl')
            bc_btag_SF_lookup = tt_cor.GetBtagSFLookup(wp='M',year=year, method='deepJet_comb')

            raw_met = met
            cleanedJets['pt_orig'] = cleanedJets.pt     # NECESSARY FOR MET CORRECTIONS LATER 
            rho_jagged = ak.ones_like(cleanedJets.pt) * events.fixedGridRhoFastjetAll
            cleanedJets = ak.with_field(cleanedJets, rho_jagged, "Rho")
            cleanedJets["pt_gen"] = ak.values_astype(ak.fill_none(cleanedJets.matched_gen.pt, 0), np.float32)
            cleanedJets = tt_cor.ApplyJetCorrections(year, corr_type='jets', isData=isData, era=run_era).build(cleanedJets)
            corrected_met = tt_cor.ApplyJetCorrections(year, corr_type='met', isData=isData, era=run_era).build(MET=raw_met, corrected_jets=cleanedJets)

        ######### The rest of the processor is inside this loop over systs that affect object kinematics  ###########
        for kinematic_var in kinematic_variations: 

            if isData: 
                cleanedJets['isGood'] = tt_os.is_pres_jet(cleanedJets)
                goodJets =  cleanedJets[cleanedJets.isGood]

                jet_veto_map = tt_cor.ApplyJetVetoMaps(goodJets, year)

                njets = ak.num(goodJets)
                jets_sorted = goodJets[ak.argsort(goodJets.pt, axis=-1,ascending=False)]
                jets_sorted = ak.pad_none(jets_sorted, 1)
                j0 = jets_sorted[:,0]
                isBtagJetsMedium = (goodJets.btagDeepFlavB > btagwpm)
                nbtagsm = ak.num(goodJets[isBtagJetsMedium])


            elif not isData:

                weights_obj_base_for_kinematic_syst = copy.deepcopy(weights_obj_base)

                if kinematic_var == 'nominal': 
                    cleanedJets['isGood'] = tt_os.is_pres_jet(cleanedJets)
                    goodJets =  cleanedJets[cleanedJets.isGood]
                    met = corrected_met

                elif 'METunclust' in kinematic_var: #  (kinematic_var == 'METunclustUp') or (kinematic_var == 'METunclustDown'):
                    print(f"\n doing METunclust systmatic")
                    cleanedJets['isGood'] = tt_os.is_pres_jet(cleanedJets)
                    goodJets =  cleanedJets[cleanedJets.isGood]
                    met = tt_cor.GetMETunclust(events, year, original_obj=corrected_met, syst_var=kinematic_var)

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

                        weights_obj_base_for_kinematic_syst.add(f'btagSF{b_syst}', events.nom, event_weight_up/btag_eventweight, event_weight_down/btag_eventweight)


            ######## Create objects for dense axes ##########
            bjets = goodJets[isBtagJetsMedium] 
            bjets_sorted = bjets[ak.argsort(bjets.pt, axis=-1,ascending=False)] 
            bjets_padded = ak.pad_none(bjets_sorted, 2)

            b0 = ak.fill_none(bjets_padded[:, 0], 0)
            b1 = ak.fill_none(bjets_padded[:, 1], 0)

            def get_sum_mass(obj_list):
                # Manually calculate px, py, pz, and energy for each object
                total_px = 0
                total_py = 0
                total_pz = 0
                total_e  = 0
                
                for obj in obj_list:
                    # Calculate components from pt, eta, phi, mass
                    px = obj.pt * np.cos(obj.phi)
                    py = obj.pt * np.sin(obj.phi)
                    pz = obj.pt * np.sinh(obj.eta)
                    # Energy = sqrt(p^2 + m^2)
                    p2 = px**2 + py**2 + pz**2
                    energy = np.sqrt(p2 + obj.mass**2)
                    
                    total_px = total_px + px
                    total_py = total_py + py
                    total_pz = total_pz + pz
                    total_e  = total_e + energy
                
                # Invariant mass: M = sqrt(E^2 - p^2)
                m2 = total_e**2 - (total_px**2 + total_py**2 + total_pz**2)
                return np.sqrt(np.maximum(0, m2))

            # Use the function for your 4-object system
            mllbb = get_sum_mass([ln, lp, b0, b1])
            mll   = get_sum_mass([ln, lp])
            mbb   = get_sum_mass([b0, b1])

            ######## HEM veto ########
            HEM_veto_mask, HEM_event_weight = tt_es.getHemMask(events, year, isData, jets=goodJets)
            if not isData:
                weights_obj_base_for_kinematic_syst.add(f'HEM', HEM_event_weight)
            
            ######### Store boolean masks with PackedSelection ##########
            selections = PackedSelection(dtype='uint64')

            if isData: 
                selections.add('is_good_lumi', lumi_mask)

            selections.add('pass_trg', pass_trg)

            if isData: 
                selections.add("ee",  (events.is_ee & events.is2los & pass_trg))    # data always has trig pass requirement
            else: 
                selections.add("ee",  (events.is_ee & events.is2los))               # MC for ee has efficiencies, so no pass_trg requirement
                
            selections.add("em",  (events.is_em & events.is2los & pass_trg))        # MC for emu has SF, so use pass_trg requirement 
            selections.add("mm",  (events.is_mm & events.is2los & pass_trg))        # MC for mumu has SF, so use pass_trg requirement

            selections.add("DYmask", (events.mllDYmask))
            selections.add("minMETpt", (met.pt > 60))

            selections.add('jetvetomap', (jet_veto_map == 0))
            selections.add('HEMvetomap', HEM_veto_mask)

            selections.add('bmask_exactly0med', (nbtagsm==0))
            selections.add('bmask_exactly1med', (nbtagsm==1))
            selections.add('bmask_exactly2med', (nbtagsm==2))
            selections.add('bmask_atleast2med', (nbtagsm>=2))

            selections.add("exactly_0j", (njets==0))
            selections.add("exactly_1j", (njets==1))
            selections.add("exactly_2j", (njets==2))
            selections.add("exactly_3j", (njets==3))

            selections.add("atleast_1j", (njets>=1))
            selections.add("atleast_2j", (njets>=2))
            selections.add("atleast_3j", (njets>=3))
            selections.add("atleast_4j", (njets>=4))

            ######### Fill dense axes variables ##########
            features = torch.from_numpy(np.concatenate([
                [ln.pt.to_numpy()],
                [lp.pt.to_numpy()],
                [l0.pt.to_numpy()],
                [l1.pt.to_numpy()],
                [ln.eta.to_numpy()],
                [lp.eta.to_numpy()],
                [l0.eta.to_numpy()],
                [l1.eta.to_numpy()],
                [(ln+lp).pt.to_numpy()],
                [mll.to_numpy()], 
                [ln.delta_phi(lp).to_numpy()],
                [abs(ln.eta - lp.eta).to_numpy()],
                [b0.pt.to_numpy()],
                [b1.pt.to_numpy()],
                [b0.eta.to_numpy()],
                [b1.eta.to_numpy()],
                [(b0+b1).pt.to_numpy()],
                [mbb.to_numpy()],
                # [b0.delta_phi(b1).to_numpy()],
                # [abs(b0.eta - b1.eta).to_numpy()],
                # [b0.delta_phi(lp).to_numpy()],
                # [abs(b0.eta - lp.eta).to_numpy()],
                # [b0.delta_phi(ln).to_numpy()],
                # [abs(b0.eta - ln.eta).to_numpy()],
                # [b1.delta_phi(lp).to_numpy()],
                # [abs(b1.eta - lp.eta).to_numpy()],
                # [b1.delta_phi(ln).to_numpy()],
                # [abs(b1.eta - ln.eta).to_numpy()],
                [mllbb.to_numpy()],
                [njets.to_numpy()],
                [np.concatenate([[(b0+b1).pt.to_numpy()],
                                 [(l0+l1).pt.to_numpy()],
                                 [(l0 + b0).pt.to_numpy()]],
                                 axis=0).max(0)],
            ]).astype(self._dtype).T)

        fit_coefs = torch.from_numpy(eft_coeffs * norm)

        training, validation = torch.utils.data.random_split(torch.utils.data.TensorDataset(features, fit_coefs), [0.8, 0.2], generator=torch.Generator().manual_seed(42))

        torch.save(torch.utils.data.TensorDataset(training[:][0], training[:][1]), f'{self._outname}/training/{int(time.time())}{random.randint(1000000,9999999)}.p')
        torch.save(torch.utils.data.TensorDataset(validation[:][0], validation[:][1]), f'{self._outname}/validation/{int(time.time())}{random.randint(1000000,9999999)}.p')

        output = {'out':  out}

        return output
        
    def postprocess(self, accumulator):
        return accumulator

