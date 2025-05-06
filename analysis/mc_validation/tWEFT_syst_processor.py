#!/usr/bin/env python
import coffea
import numpy as np
import awkward as ak
import json
import hist

from mt2 import mt2
from mt2 import mt2_arxiv

from coffea import processor
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods import vector

# silence warnings due to using NanoGEN instead of full NanoAOD
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
NanoAODSchema.warn_missing_crossrefs = False

from topcoffea.modules.histEFT import HistEFT
import topcoffea.modules.eft_helper as efth
import topcoffea.modules.corrections as tc_cor

from ttbarEFT.modules.analysis_tools import make_mt2, get_lumi
from ttbarEFT.modules.axes import tW_axes as axes_info
import ttbarEFT.modules.selections as selec

np.seterr(divide='ignore', invalid='ignore', over='ignore')

class AnalysisProcessor(processor.ProcessorABC):
    def __init__(self, samples, wc_names_lst=[], hist_lst = None, dtype=np.float32, do_errors=False, do_systs=False):
        self._samples = samples
        self._wc_names_lst = wc_names_lst

        # set the booleans
        self._dtype = dtype
        self._do_errors = do_errors # Whether to calculate and store the w**2 coefficients
        self._do_systematics = do_systs # Whether to process systematic samples

        # create axes
        proc_axis = hist.axis.StrCategory([], name="process", growth=True)
        chan_axis = hist.axis.StrCategory([], name="channel", growth=True)
        syst_axis = hist.axis.StrCategory([], name="systematic", label=r"Systematic Uncertainty", growth=True)

        # fill histograms using info from axes.py 
        histograms = {}
        for name, info in axes_info.items():
            if "variable" in info: 
                dense_axis = hist.axis.Variable(info["variable"], name=name, label=info["label"])
                sumw2_axis = hist.axis.Variable(info["variable"], name=name+"_sumw2", label=info["label"]+" sum of w^2")
            else:
                dense_axis = hist.axis.Regular(*info["regular"], name=name, label=info["label"])
                sum2w_axis = hist.axis.Regular(*info["regular"], name=name+"_sumw2", label=info["label"] + " sum of w^2")

            histograms[name] = HistEFT(
                proc_axis, 
                syst_axis,
                dense_axis,
                wc_names = wc_names_lst, 
                label=r"Events",
            )

            histograms[name+"_sumw2"] = HistEFT(
                proc_axis, 
                syst_axis,
                sum2w_axis,
                wc_names = wc_names_lst, 
                label=r"Events",
            )

        self._accumulator = histograms

        # set the list of hists to fill
        if hist_lst is None:
            self._hist_lst = list(self._accumulator.keys()) #fill all hists if not specified
        else:
            for hist_to_include in hist_lst:
                if hist_to_include not in self._accumulator.keys():
                    raise Exception(f"Error: Cannot specify hist \"{hist_to_include}\", it is not defined in the processor.")
            self._hist_lst = hist_lst 

        # print out basic info before running over files
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
        isData          = self._samples[dataset]["isData"]
        hist_axis_name  = self._samples[dataset]["histAxisName"]
        year            = self._samples[dataset]['year']
        xsec            = self._samples[dataset]['xsec']
        sow             = self._samples[dataset]['nSumOfWeights']

        isEFT = hasattr(events, "EFTfitCoefficients")

        assert not (isEFT and isData), f"isEFT and isData cannot both be True. Check input samples."

        # Get up down weights from input dict
        if (self._do_systematics and not isData):
            sow_ISRUp          = self._samples[dataset]["nSumOfWeights"]
            sow_ISRDown        = self._samples[dataset]["nSumOfWeights"]
            sow_FSRUp          = self._samples[dataset]["nSumOfWeights"]
            sow_FSRDown        = self._samples[dataset]["nSumOfWeights"]
            sow_renormUp       = self._samples[dataset]["nSumOfWeights"]
            sow_renormDown     = self._samples[dataset]["nSumOfWeights"]
            sow_factUp         = self._samples[dataset]["nSumOfWeights"]
            sow_factDown       = self._samples[dataset]["nSumOfWeights"]
            sow_renormfactUp   = self._samples[dataset]["nSumOfWeights"]
            sow_renormfactDown = self._samples[dataset]["nSumOfWeights"]

            # sow_ISRUp          = self._samples[dataset]["nSumOfWeights_ISRUp"]
            # sow_ISRDown        = self._samples[dataset]["nSumOfWeights_ISRDown"]
            # sow_FSRUp          = self._samples[dataset]["nSumOfWeights_FSRUp"]
            # sow_FSRDown        = self._samples[dataset]["nSumOfWeights_FSRDown"]
            # sow_renormUp       = self._samples[dataset]["nSumOfWeights_renormUp"]
            # sow_renormDown     = self._samples[dataset]["nSumOfWeights_renormDown"]
            # sow_factUp         = self._samples[dataset]["nSumOfWeights_factUp"]
            # sow_factDown       = self._samples[dataset]["nSumOfWeights_factDown"]
            # sow_renormfactUp   = self._samples[dataset]["nSumOfWeights_renormfactUp"]
            # sow_renormfactDown = self._samples[dataset]["nSumOfWeights_renormfactDown"]

        else:
            sow_ISRUp          = -1
            sow_ISRDown        = -1
            sow_FSRUp          = -1
            sow_FSRDown        = -1
            sow_renormUp       = -1
            sow_renormDown     = -1
            sow_factUp         = -1
            sow_factDown       = -1
            sow_renormfactUp   = -1
            sow_renormfactDown = -1


        ########## Initialize Objects ##########
        genpart = events.GenPart
        is_final_mask = genpart.hasFlags(["fromHardProcess","isLastCopy"])

        ele  = genpart[is_final_mask & (abs(genpart.pdgId) == 11)]
        mu   = genpart[is_final_mask & (abs(genpart.pdgId) == 13)]
        nu_ele = genpart[is_final_mask & (abs(genpart.pdgId) == 12)]
        nu_mu = genpart[is_final_mask & (abs(genpart.pdgId) == 14)]
        nu = ak.concatenate([nu_ele,nu_mu],axis=1)

        jets = events.GenJet
        gen_top = genpart[is_final_mask & (abs(genpart.pdgId) == 6)]
        met = events.GenMET

        # An array of length events that is just 1 for each event
        events.nom = ak.ones_like(events.GenMET)


        ########## EFT Coefficients ##########
        eft_coeffs = ak.to_numpy(events["EFTfitCoefficients"]) if hasattr(events, "EFTfitCoefficients") else None
        # Check to see if the ordering of WCs for this sample matches what want
        if eft_coeffs is not None:
            if self._samples[dataset]["WCnames"] != self._wc_names_lst:
                eft_coeffs = efth.remap_coeffs(self._samples[dataset]["WCnames"], self._wc_names_lst, eft_coeffs)
        eft_w2_coeffs = efth.calc_w2_coeffs(eft_coeffs,self._dtype) if (self._do_errors and eft_coeffs is not None) else None


        ########## Object Selection ##########
        # currently, these are returned sorted by pt, might need to later change to sort in the processor
        leps = selec.gen_lepton_selections(ele, mu)
        jets = selec.gen_jets_selection(jets, leps, nu)

        nleps = ak.num(leps)
        ntops = ak.num(gen_top)
        # l0 = leps[ak.argmax(leps.pt, axis=-1, keepdims=True)]
        leps_padded = ak.pad_none(leps, 2)
        l0 = leps_padded[:,0]
        l1 = leps_padded[:,1]

        ########## Systematics ##########

        obj_correction_syst_lst = []
        data_syst_lst = []

        # wgt_correction_syst_lst = [
        #     "lepSF_muonUp","lepSF_muonDown","lepSF_elecUp","lepSF_elecDown",f"btagSFbc_{year}Up",f"btagSFbc_{year}Down","btagSFbc_corrUp","btagSFbc_corrDown",f"btagSFlight_{year}Up",f"btagSFlight_{year}Down","btagSFlight_corrUp","btagSFlight_corrDown","PUUp","PUDown","PreFiringUp","PreFiringDown",f"triggerSF_{year}Up",f"triggerSF_{year}Down", # Exp systs
        #     "FSRUp","FSRDown","ISRUp","ISRDown","renormUp","renormDown","factUp","factDown", # Theory systs
        #     ]

        # Define the lists of systematics we include
        wgt_correction_syst_lst = [
            "FSRUp","FSRDown","ISRUp","ISRDown","renormUp","renormDown","factUp","factDown", "renormfactDown", "renormfactUp",# Theory systs
            ]

        # We only calculate these values if not isData
        weights_obj_base = coffea.analysis_tools.Weights(len(events),storeIndividual=True)
        if not isData: 
            if eft_coeffs is None:
                genw = events["genWeight"]
            else:
                genw = np.ones_like(events['event'])

            # Normalize by (xsec/sow)*genw where genw is 1 for EFT samples
            # Note that for theory systs, will need to multiply by sow/sow_wgtUP to get (xsec/sow_wgtUp)*genw and same for Down
            lumi = 1000.0*get_lumi(year)
            weights_obj_base.add("norm",(xsec/sow)*genw*lumi)

            # Attach PS weights (ISR/FSR) and scale weights (renormalization/factorization) and PDF weights
            tc_cor.AttachPSWeights(events)
            tc_cor.AttachScaleWeights(events)

            # FSR/ISR weights -- corrections come from AttachPSWeights
            weights_obj_base.add('ISR', events.nom, events.ISRUp*(sow/sow_ISRUp), events.ISRDown*(sow/sow_ISRDown))
            weights_obj_base.add('FSR', events.nom, events.FSRUp*(sow/sow_FSRUp), events.FSRDown*(sow/sow_FSRDown))
            # renorm/fact scale  -- corrections come from AttachScaleWeights
            weights_obj_base.add('renorm', events.nom, events.renormUp*(sow/sow_renormUp), events.renormDown*(sow/sow_renormDown))
            weights_obj_base.add('fact', events.nom, events.factUp*(sow/sow_factUp), events.factDown*(sow/sow_factDown))
            weights_obj_base.add('renormfact', events.nom, events.renormfactUp*(sow/sow_renormfactUp), events.renormfactDown*(sow/sow_renormfactDown))


        ########## The rest of the processor is inside this loop over systs that affect object kinematics ##########
        
        # If we're doing systematics and this isn't data, we will loop over the obj_correction_syst_lst list
        if self._do_systematics and not isData: 
            syst_var_list = ["nominal"] + obj_correction_syst_lst
        else: 
            syst_var_list = ["nominal"]

        # Loop over the list of systematic variations we've constructed
        for syst_var in syst_var_list:
            # Make a copy of the base weights object, so that each time through the loop we do not double count systs
            # In this loop over systs that impact kinematics, we will add to the weights objects the SFs that depend on the object kinematics
            weights_obj_base_for_kinematic_syst = copy.deepcopy(weights_obj_base)

            #### Stuff that effects jets will go here once we include regular systematics #### 
            # that will result in a new jets object that has corrections applied
            # that object will then be used below in `njets`, `j0`, `ht`, etc. 

            njets = ak.num(jets)
            j0 = jets[ak.argmax(jets.pt, axis=-1, keepdims=True)]

            #LATER: ########## Add variables into event object so that they persist ##########
            #LATER: ########## Event weights that do not depend on the lep cat ##########
            #LATER: ########## Event weights that do depend on the lep cat ##########

            weights_dict={} #in the TOP-22-006 processor, the weights_dict has a different entry for each channel, and the following section is in a loop over channels

            weights_dict['tW'].copy.deepcopy(weights_obj_base_for_kinematic_syst)


            ########## Event Selection Masks ##########
            at_least_two_leps = ak.fill_none(nleps>=2,False)
            at_least_one_jet = ak.fill_none(njets>=1,False)


            ########## Store boolean masks with PackedSelection ##########
            selections = PackedSelection()
            selections.add('2l', at_least_two_leps)
            selections.add('1j', at_least_one_jet)

            event_selection_mask = selections.all('2l', '1j')


            ######### Variables for the dense axes of the hists ##########
            ht = ak.sum(jets.pt, axis=-1)
            mt2_var = make_mt2(l0, l1, met)

            dr_leps = l0.delta_r(l1)
            mll = (l0+l1).mass

            counts = np.ones_like(events['event'])

            # Variables we will loop over when filling hists
            varnames = {}
            varnames["njets"] = njets
            varnames["nleps"] = nleps 
            varnames["ntops"] = ntops
            varnames["top_pt"] = gen_top.pt
            varnames["l0pt"] = ak.flatten(l0.pt) 
            varnames["dr_leps"] = dr_leps
            varnames["mll"] = mll
            varnames["mt2"] = mt2_var
            varnames["sow"] = counts

            ########## Fill Histos ##########

            for dense_axis_name, dense_axis_vals in varnmaes.items():
                if dense_axis_name not in self._hist_lst: continue

                wgt_var_lst = ["nominal"]
                if self._do_systematics:
                    if not isData: 
                        if (syst_var != "nominal"):
                            # In this case, we are dealing with systs that change the kinematics of the objs (e.g. JES)
                            # So we don't want to loop over up/down weight variations here
                            wgt_var_lst = [syst_var]
                        else:
                            # Otherwise we want to loop over the up/down weight variations
                            wgt_var_lst = wgt_var_lst + wgt_correction_syst_lst

                    else:
                        # This is data, so we want to loop over just up/down variations relevant for data (i.e. FF up and down)
                        wgt_var_lst = wgt_var_lst + data_syst_lst

                # Loop over the systematics
                for wgt_fluct in wgt_var_lst:
                    weights_object = weights_dict["tW"] # usually this is inside a loop over nlep categories
                    if (wgt_fluct == "nominal") or (wgt_fluct in obj_correction_syst_lst):
                        # In the case of "nominal", or the jet energy systematics, no weight systematic variation is used
                        weight = weights_object.weight(None)
                    else: 
                        # Otherwise get the weight from the Weights object
                        if wgt_fluct in weights_object.variations:
                            weight = weights_object.weight(wgt_fluct)
                        else:
                            # Note in this case there is no up/down fluct for this cateogry, so we don't want to fill a hist for it
                            continue

                    weights_flat = weight[all_cuts_mask]
                    eft_coeffs_cut = eft_coeffs[all_cuts_mask] if eft_coeffs is not None else None

                    # Fill the histos
                    axes_fill_info_dict = {
                        dense_axis_name : dense_axis_vals[event_selection_mask],
                        "process"       : histAxisName,
                        "systematic"    : wgt_fluct,
                        "weight"        : weights_flat,
                        "eft_coeff"     : eft_coeffs_cut,
                    }

                    hout[dense_axis_name].fill(**axes_fill_info_dict)


                    axes_sumw2_fill_info_dict = {
                        dense_axis_name+"_sumw2" : dense_axis_vals[all_cuts_mask],
                        "process"       : histAxisName,
                        "systematic"    : wgt_fluct,
                        "weight"        : np.square(weights_flat),
                        "eft_coeff"     : eft_coeffs_cut,
                    }
                    hout[dense_axis_name+"_sumw2"].fill(**axes_sumw2_fill_info_dict)

            hout = self.accumulator

    def postprocess(self, accumulator):
        return accumulator