import pickle
import gzip
import numpy as np
import boost_histogram as bh
import uproot
import hist
import os
import re
import json
import yaml
import time

from collections import defaultdict
from topcoffea.modules.utils import regex_match
from ttbarEFT.modules.paths import ttbarEFT_path

PRECISION = 6   # Decimal point precision in the text datacard output


def to_hist(arr,name,zero_wgts=False):
    """
        Converts a numpy array into a hist.Hist object suitable for being written to a root file by
        uproot. 
        arr[0]: bin contents (sum of weights)
        arr[1]: stat uncertaitny
        zero_wgts=True means the resulting histogram will be created with bin errors set to 0

    # NOTE:
    #   If we don't instantiate a new np.array here, then clipped will store a reference to the
    #   sub-array arr and when we modify clipped, it will propagate back to arr as well!
    """
    clipped = []
    for i in range(2):  # first entry is sum(weight), second entry is sum(weight^2)
        if arr[i] is not None:
            if "PDF" in name:   # don't remove the 0th bin from PDF, since it's calculated separately it's already removed 
                clipped.append(np.array(arr[i]))
            else: 
                clipped.append(np.array(arr[i][1:]))  # Strip off the underoverflow bin
            # clipped.append(np.array(arr[i][1:]))  # Strip off the underoverflow bin
        else:
            clipped[i] = None

    nbins = len(clipped[0])
    h = hist.Hist(hist.axis.Regular(nbins,0,nbins,name=name),storage=bh.storage.Weight())   #bh.storage.Weight() stores the value and variance of each bin
    if zero_wgts:
        h[...] = np.stack([clipped[0],np.zeros_like(clipped[0])],axis=-1) # Set the bin errors all to 0
    else:
        h[...] = np.stack([clipped[0], clipped[1]],axis=-1)
    return h


def add_sumw2_stub(eval_d, sumw2=None):
    # compatibility with the coffea.Hist. add a row of zeros for sumw2.
    eval_d2 = {}
    for k, v in eval_d.items():
        if sumw2 is None:
            eval_d2[k] = np.stack((v, np.broadcast_to(np.zeros((1,)), len(v))))    # If no sumw2, just stack with zeros
        else:
            if k in sumw2:
                variance = sumw2[k]
            else:   #Find the key in sumw2 that matches everything EXCEPT systematic
                try:
                    # Look for the entry in sumw2 that matches the process/channel
                    var_key = next(sk for sk in sumw2.keys() if all(getattr(sk, a) == getattr(k, a) for a in sk._fields))
                    variance = sumw2[var_key]
                except (StopIteration, AttributeError):
                    variance = np.zeros_like(v)
            
            eval_d2[k] = np.stack((v, variance))
    
    return eval_d2


class RateSystematic():
    def __init__(self,name,**kwargs):
        self.all = kwargs.pop("all",False)      # If true, this syst applies to all processes
        if self.all:
            try:
                self.all_unc = kwargs.pop("unc")
            except KeyError:
                msg = "Missing 'unc' argument. Must specify an uncertainty when using the 'all' option"
                raise KeyError(msg)
        self.name = name
        self.corrs = {}  # keys are the name of processes and values are the corresponding unc.

    def has_process(self,p):
        return self.all or (p in self.corrs)

    def add_process(self,p,v=None):
        if self.all:
            raise KeyError("Can't add a correlated process for systematic defined with the 'all' option")
        self.corrs[p] = v

    # TODO: This needs to be given a better name
    # Returns the corresponding unc. (i.e. kappa values) that have been associated with a particular process
    # Note: The return value should be as a string
    def get_process(self,p):
        if self.all:
            return self.all_unc
        if self.has_process(p):
            return self.corrs[p]
        else:
            # This is the case for a systematic that doesn't apply to the specified process
            return '-'


class DatacardMaker():
    # Note:
    #   Care must be taken with regards to the underscores, due to 'nonprompt', 'data', and 'flips'
    GROUP = {
        'tt': [
            'TT01j2lmtt0to700_',
            'TT01j2lmtt700to900_',
            'TT01j2lmtt900toInf_'
            ],
        'tW': [
            'TW_NoFullyHadronicDecays_',
            ],
        'DY': [
            'DY10to50_',
            'DY50_',
            ],
        'Others':[
            'WJetsToLNu_',     
            'WWTo2L2Nu_',
            'WZTo3LNu_',
            'ZZTo4L_',
            'TTGJets_',
            'ttW_',
            'ttZ_',
            ],
        'data_obs': ['data'],
    }

    # Controls how we rebin the dense axis of the corresponding distribution
    # TODO update this to read the yaml correctly
    # BINNING = {}
    # for name, value in axes_info.items():
    #     if "variable" in value:
    #         BINNING[name] = value["variable"]

    YEARS = ["UL16","UL16APV","UL17","UL18"]
    SYST_YEARS = ["2016","2016APV","2017","2018"]

    # TODO update this to my namign 
    FNAME_TEMPLATE = "TESTING-ttbar-{cat}_{kmvar}.{ext}"
    # SIGNALS = set(["ttH","tllq","ttll","ttlnu","tHq","tttt"])
    SIGNALS = set(["tt"])

    @classmethod
    def get_year(cls,s):
        """
            Attempt to return the year of the process or systematic string
        """
        for yr in cls.YEARS:
            if s.endswith(yr): return yr
        for yr in cls.SYST_YEARS:
            if s.endswith(yr+"Up"): return yr
            if s.endswith(yr+"Down"): return yr
        return None

    @classmethod
    def strip_fluctuation(cls,s):
        if s.endswith("Down"):
            return s[:-4]
        elif s.endswith("Up"):
            return s[:-2]
        else: 
            return s

    @classmethod
    def strip_year(cls,s):
        for yr in cls.YEARS:
            s = s.replace(yr,"")
        for yr in cls.SYST_YEARS:
            s = s.replace(f"_{yr}","")  # Note the underscore
        return s

    @classmethod
    def is_signal(cls,s):
        s = cls.get_process(s)
        return (s in cls.SIGNALS)

    @classmethod
    def is_per_year_systematic(cls,s):
        end_chks = [
            "_2016APVUp","_2016Up","_2017Up","_2018Up",
            "_2016APVDown","_2016Down","_2017Down","_2018Down",
        ]
        return any([s.endswith(x) for x in end_chks])

    @classmethod
    def is_eft_term(cls,s):
        """ Check if string corresponds an EFT process term after decomposition."""
        chks = ["_lin_","_quad_"]
        return any([x in s for x in chks])

    @classmethod
    def get_process(cls,s):
        """ Strips off the year designation from a process name, can also be used for decomposed terms."""
        for yr in cls.YEARS:
            if s.endswith(yr):
                s = s.replace(yr,"")
        if cls.is_eft_term(s):
            # For now we can simply split on first underscore, if signal process names get underscores
            #   will need to update this to be smarter
            s = s.split("_",1)[0]
        if "_" in s:
            s = s.rsplit("_",1)[0]
        return s

    # works for ttbarEFT
    @classmethod
    def get_jet_mults(cls,s):
        """
            Returns the njet and bjet multiplicities based on the string passed to it in (j,b) order.
            For the regular expression, group 1 matches 'njet_bjet', group 2 matches 'bjet_njet'
            group 3 matches '_njet'.
        """
        # rgx = re.compile(r"(_[2-7]j_[1-2]b)|(_[1-2]b_[1-7]j)|(_[1-7]j$)")
        rgx = re.compile(r"(_[1-7]j_[1-2]b)|(_[1-2]b_[1-7]j)|(_[1-7]j$)")

        m = rgx.search(s.replace('_fwd', ''))
        if m.group(1) and m.group(2) is None and m.group(3) is None:
            # The order is '_Nj_Mb'
            _,j,b = m.group(1).split("_")
        elif m.group(1) is None and m.group(2) and m.group(3) is None:
            # The order is '_Nb_Mj'
            _,b,j = m.group(2).split("_")
        elif m.group(1) is None and m.group(2) is None and m.group(3):
            # This occurs when the string ends in '_Mj' and doesn't have a bjet multiplicity
            b = None
            j = m.group(3).replace("_","")
        else:
            raise ValueError(f"Unable to find rgx match in string {s}")
        j = int(j.replace("j",""))
        if b is not None:
            b = int(b.replace("b",""))
        return (j,b)

    @classmethod
    def get_lep_mult(cls,s):
        """ Returns the lepton multiplicity based on the string passed to it."""
        return 2

    @classmethod
    def get_processes_by_years(cls,h):
        """
            Reads the 'process' sparse axis of a histogram and returns a dictionary that maps stripped
            process names to the list of sparse axis categories it came from.
        """
        r = defaultdict(lambda: [])
        for x in h.axes["process"]:
            p = cls.get_process(x)
            r[p].append(cls.get_process(x))
        return r

    def __init__(self,pkl_path,**kwargs):
        self.year_lst        = kwargs.pop("year_lst",[])
        self.do_sm           = kwargs.pop("do_sm",False)
        self.do_nuisance     = kwargs.pop("do_nuisance",False)
        self.drop_syst       = kwargs.pop("drop_syst",[])
        self.out_dir         = kwargs.pop("out_dir",".")
        self.var_lst         = kwargs.pop("var_lst",[])
        self.do_mc_stat      = kwargs.pop("do_mc_stat",False)
        self.coeffs          = kwargs.pop("wcs",[])
        self.use_real_data   = kwargs.pop("unblind",False)
        self.verbose         = kwargs.pop("verbose",True)
        self.use_AAC         = kwargs.pop("use_AAC",False)
        self.wc_scalings     = kwargs.pop("wc_scalings",[])
        self.scalings        = []

        # we don't need the rotations for ttbar SMFEFTsim
        self.rotate = {}
        
        ### Initial Settings ###
        if self.year_lst:
            for yr in self.year_lst:
                if not yr in self.YEARS:
                    raise ValueError(f"Invalid year choice '{yr}', should be empty if running over all years or one of: {self.YEARS}")

        # get wc ranges from json
        with open(ttbarEFT_path("params/wc_ranges.json"), "r") as wc_ranges_json:
            self.wc_ranges = json.load(wc_ranges_json)

        rate_syst_path = kwargs.pop("rate_systs_path", "params/rate_systs.json")
        self.rate_systs = self.load_systematics(rate_syst_path)

        # Samples to be excluded from the datacard, should correspond to names before group_processes is run
        self.ignore = [
            # "DYJetsToLL", "DY10to50", "DY50",
            # data,
        ]

        # Since we're just going to generate Asimov data, this lets us drop the real data histograms for a minor speed-up
        if not self.use_real_data:
            print(f"dropping real data histograms since we're making Asimov")
            self.ignore.append("data")
        extra_ignore = kwargs.pop("ignore",[])
        if extra_ignore:
            print(f"Adding processes to ignore: {extra_ignore}")
        self.ignore.extend(extra_ignore)

        # the processes in this dictionary should match the naming AFTER grouping
        self.syst_skip = {
            "PDF": ["tW", "DY", "Others"],
            "hdamp": ["tW", "DY", "Others"],
        }

        """
        For now, we leave this as a hardcoded thing, a bit tedious but it works
        Note: If not explicitly listed, it is assumed that all years should be uncorrelated
        Note: It is important to list the correlations for ALL years in which it is relevant, so
              for example, if a systematic is correlated in 2016, 2016APV, and 2017, there needs
              to be an entry for all three years and the list corresponding to each entry needs to
              be consistent (i.e. contain the other correlated years) across all three entries
        Note: As a final note, the actual systematic that appears in the datacards will be just one
              of the set to be correlated. So for example, if a systematic is correlated over 2016
              and 2016APV, then either the 2016 or 2016APV version will appear in the datacard, but
              not both. Typically, the one that remains will be the 2016 version as that's the one
              that gets handled first in the loop, but it would be different if we processed things
              in a different order
        """
        self.syst_year_corr = {
            # Example of correlated for only 2016 and 2016APV
            # "FFcloseEl": {"2016": ["2016APV"], "2016APV": ["2016"]},

            # Example of correlated over 2016, 2016APV, 2017, and 2018
            # "FFcloseEl": {"2016": ["2016APV","2017","2018"], "2016APV": ["2016","2017","2018"], "2017": ["2016","2016APV","2018"], "2018": ["2016","2016APV","2017"]},
        }

        """
        Defines which systematics should be decorrelated in the self.analysis() step. 
            Each key should match (exactly) a particular systematic. 
            The list for each systematic specifies which processes should remain remain correlated or not.
        Note: A given process should appear AT MOST once in the "matches" list for a given systematic
            grouping. If a process has an associated systematic, but doesn't match any of the
            groups, then it will retain its original systematic name (i.e. all unmatched
            processes will remain correlated).
        Note: For the special case where the group name is an empty string the systematic will
            instead have the matched process' name appended to it, meaning that all matched
            processes will be decorrelated!
        Note: Since the decorrelation happens during the self.analysis() step, the matched names
            should correspond to the renamed/re-grouped processes, e.g. use "Diboson" instead of
            "ZZ","WZ","WW". 
        """
        self.syst_shape_decorrelate = {
            "renorm": [{
                "matches": ["tt", "tW", "DY", "Others"],
                "group": "",
            }],
            "fact": [{
                "matches": ["tt", "tW", "DY", "Others"],
                "group": "",
            }]
        }
        
        if extra_ignore:
            print(f"Adding processes to ignore: {extra_ignore}")
        self.ignore.extend(extra_ignore)

        self.tolerance = 1e-4
        self.hists = None

        tic = time.time()
        self.read(pkl_path)
        dt = time.time() - tic
        print(f"Total Read+Prune Time: {dt:.2f} s")
        print (f"Saving output to {os.path.realpath(self.out_dir)}")
        
        
    def load_systematics(self,rs_fpath):
        """
            Parse out the correlated and decorrelated systematics from rate_systs.json and
            missing_parton.root files.
        """
        rate_systs = {}
        if not self.do_nuisance:
            return rate_systs
        fpath = ttbarEFT_path(rs_fpath)
        print(f"Opening: {fpath}")
        with open(fpath) as f:
            rates_json = json.load(f)
        for k1,v1 in rates_json["rate_uncertainties"].items():
            # k1 will be the name of a rate systematic, like 'lumi' or 'pdf_scale'
            if isinstance(v1,dict):
                # This is a correlated rate systematic
                syst_name = f"{k1}"
                new_syst = RateSystematic(syst_name)
                for k2,v2 in v1.items():
                    # k2 will be the name of process, like 'charge_flips' or 'ttH' or 'Diboson'
                    p = self.get_process(k2)
                    new_syst.add_process(p,v2)
                rate_systs[k1] = new_syst
            else:
                # The systematic gets applied to everything
                syst_name = f"{k1}"
                new_syst = RateSystematic(syst_name,all=True,unc=v1)
                rate_systs[k1] = new_syst

        # Certain rate systematics are only correlated between subsets of processes
        to_remove = set()
        for p,corr_systs in rates_json["correlations"].items():
            # 'p' is the name of a process and 'corr_systs' is a dictionary defining which rate systematic
            #  needs to be decorrelated into a specific sub-group, e.g. for ttH: pdf_scale -> pdf_scale_gg
            for syst,grp in corr_systs.items():
                # 'syst' should be the name of a systematic already defined in the 'rate_systs' dictionary
                #   and 'grp' is the string we are going to differentiate it from the other variants of 'syst'
                to_remove.add(syst)
                syst_name = f"{syst}_{grp}"
                if not syst_name in rate_systs:
                    rate_systs[syst_name] = RateSystematic(syst_name)
                if not rate_systs[syst].has_process(p):
                    print(f"Warning: No process {p} found for {syst} systematic")
                    continue
                unc = rate_systs[syst].get_process(p)
                rate_systs[syst_name].add_process(p,unc)
        # Now lets remove the original systematics which we decorrelated into sub-groups
        for syst in to_remove:
            rate_systs.pop(syst)

        return rate_systs
    
    def read(self,fpath):
        """
            Input should be a file path to a pkl file containing histograms produced by the topeft.py
            processor. The histograms are extracted and then pre-processed to remove / group / scale
            various sparse axes categories.
        """
        print(f"Opening: {fpath}")
        tic = time.time()
        self.hists = pickle.load(gzip.open(fpath))
        dt = time.time() - tic
        print(f"Pkl Open Time: {dt:.2f} s")

        print(f"km_dists in read function: {self.hists.keys()}")
        for km_dist, h in self.hists.items():
            if h.empty():
                continue
            if self.var_lst and (km_dist not in self.var_lst): 
                continue
                
            print(f"Loading: {km_dist}")
            # Remove processes that we don't include in the datacard
            to_remove = []
            for x in h.axes["process"]:
                p = self.get_process(x)
                if p in self.ignore:
                    if self.verbose:
                        print(f"Skipping (ignored): {x}")
                    to_remove.append(x)
                    continue
                if self.year_lst:
                    yr = self.get_year(x)
                    if yr not in self.year_lst:
                        if self.verbose:
                            print(f"Skipping (year): {x}")
                        to_remove.append(x)
                        continue
            to_keep = [x for x in h.axes["process"] if x not in to_remove]
            if to_keep:
                h = h[{"process": to_keep}]
            else:
                continue

            if not self.do_nuisance:
                # Remove all shape systematics
                h = h[{"systematic": "nominal"}]

            if self.drop_syst:
                to_drop = set()
                for syst in self.drop_syst:
                    # If you pass "JEC", it drops "JECUp" and "JECDown"
                    # If you pass "JECUp", it only drops that one
                    if syst.endswith("Up") or syst.endswith("Down"):
                        to_drop.add(syst)
                    else:
                        to_drop.add(f"{syst}Up")
                        to_drop.add(f"{syst}Down")
                
                to_keep = [s for s in h.axes["systematic"] if s not in to_drop]
                h = h[{"systematic": to_keep}]
    
            if km_dist == "LHEPDFweights":
                print("LHEPDFweights is the km_dist before grouping...")
            # Remove 'central', 'private', '_4F' text from process names
            grp_map_clean = {}
            for x in h.axes["process"]:
                new_name = x.replace("private", "").replace("central", "").replace("_4F", "")
                # Note: The group method expects a list of old names for each new name
                if new_name not in grp_map_clean:
                    grp_map_clean[new_name] = []
                grp_map_clean[new_name].append(x)

            # Group based on histogram type
            if hasattr(h, "group"):         #histEFT has .group
                h = h.group("process", grp_map_clean)
            else:                           # LHEPdfWeights (regular Hist) path
                h = self.group_regular_hist(h, grp_map_clean)

            h = self.group_processes(h)
            
            if km_dist == "LHEPDFweights":
                LHEPDFweights_yr_grp_map = {}
                all_procs = list(h.axes["process"])
                group_names = self.GROUP.keys()
                for g in group_names:
                    LHEPDFweights_yr_grp_map[g]=[x for x in all_procs if x.startswith(g)]
                # print(f"LHEPDFweights_yr_grp_map:{LHEPDFweights_yr_grp_map}")

                # Group years for the regular histogram
                h =self.group_regular_hist(h, LHEPDFweights_yr_grp_map)
                self.hists[km_dist] = h
                continue

            else:
                h = self.correlate_years(h)
            
            if 'systematic' in h.axes.name:
                num_systs = len(h.axes["systematic"])
                print(f"Num. Systematics: {num_systs}")
                
            self.hists[km_dist] = h
            
    def group_processes(self,h):
        """
            Groups together certain processes from the 'process' axis. We also abuse this method to
            rename specific process categories. Both of which are determined by the GROUP static data
            member.
        """
        all_procs = set(h.axes["process"])
        grp_map = {}
        for grp_name,to_grp in self.GROUP.items():
            for yr in self.YEARS:
                new_name = f"{grp_name}{yr}"
                lst = []
                for x in to_grp:
                    old_name = f"{x}{yr}"
                    if old_name in all_procs:
                        lst.append(old_name)
                        all_procs.remove(old_name)
                # Note: Some processes only exist in certain channels (e.g. flips), so we need to
                #   skip them when they don't appear in the identifiers list
                if len(lst):
                    grp_map[new_name] = lst
        # Include back in everything that wasn't specified by the initial groupings
        for x in all_procs:
            grp_map[x] = [x]
            
        if hasattr(h, "group"):
            # Use the topcoffea SparseHist logic (safe for EFT hists)
            return h.group("process", grp_map)
        else:
            # Manual logic for regular Hist (safe for LHEPdfWeights floats)
            return self.group_regular_hist(h, grp_map)
                    
    def group_regular_hist(self, h, grp_map):
        """Handle grouping for standard hist.Hist objects."""
        # Create a new histogram with the same axes but the new process category
        new_proc_axis = hist.axis.StrCategory(grp_map.keys(), name="process", growth=True)
        other_axes = [ax for ax in h.axes if ax.name != "process"]
        h_grouped = hist.Hist(new_proc_axis, *other_axes, storage=h.storage_type())

        for new_name, existing_olds in grp_map.items():
            
            summed_h = h[{"process": existing_olds}][{"process": sum}]
            idx = h_grouped.axes["process"].index(new_name)
            h_grouped.view(flow=True)[idx] = summed_h.view(flow=True)
            
        return h_grouped

    def correlate_years(self,h):
        """
            Merges together different run years, taking care to treat year-specific systematics as
            uncorrelated from one another
        """
        if not self.do_nuisance:
            # Only sum over the years, don't mess with nuisance stuff
            grp_map = defaultdict(lambda: [])
            for x in h.axes["process"]:
                p = self.get_process(x)
                grp_map[p].append(x)
            h = h.group("process", grp_map)
            return h
        # This requires some fancy footwork to make work
        print("Correlating years")

        # Need to figure out which years are actually present in the histogram
        unique_proc_years = set()
        for x in h.axes["process"]:
            yr = self.get_year(x)
            unique_proc_years.add(yr)

        already_correlated = set()  # Keeps track of which systematics have already been correlated
        for sp_key in h.categorical_keys:
            proc = sp_key.process
            syst = sp_key.systematic
            proc_year = self.get_year(proc)
            syst_year = self.get_year(syst)
            if syst_year is None:
                # This ensures that the systematic in question is a per-year systematic
                continue
            if syst in already_correlated:
                if self.verbose:
                    print(f"Skipping {syst} as it was already correlated in a previous year")
                continue
            print(f"*** syst_base before stripping: {syst}")
            syst_base = self.strip_fluctuation(syst)
            syst_base = self.strip_year(syst_base)
            corr_keys = []
            for p_yr, s_yr in zip(self.YEARS, self.SYST_YEARS):
                if p_yr not in unique_proc_years:
                    # The histogram file was generated by running over a subset of the years or we are
                    # only making cards for a certain year
                    continue
                if p_yr == proc_year:
                    # We never add self to self
                    continue
                if syst_base in self.syst_year_corr and s_yr in self.syst_year_corr[syst_base] and syst_year in self.syst_year_corr[syst_base][s_yr]:
                    # The systematic for this year needs to be correlated
                    syst_key = syst.replace(syst_year, s_yr)
                    already_correlated.add(syst_key)
                else:
                    # The systematic for this year needs to be uncorrelated
                    syst_key = "nominal"
                proc_key = proc.replace(proc_year, p_yr)

                # Construct the sparse key
                corr_key = sp_key._asdict()
                corr_key["process"] = proc_key
                corr_key["systematic"] = syst_key
                corr_key = type(sp_key)(**corr_key)
                corr_keys.append(corr_key)

            for k in corr_keys:
                h[sp_key] += h[k]

            sp_tup = tuple(sp_key)
            if self.verbose:
                print(f"{tuple(sp_tup)} -- {' + '.join(map(lambda k: str(tuple(k)), corr_keys))}")

        # Finally sum over years, since the per-year systematics only appear in a corresponding
        #   "process year", the grouping for those systematics just adds itself with nothing from
        #   the other process years
        grp_map = defaultdict(lambda: [])
        for x in h.axes["process"]:
            p = self.get_process(x)
            grp_map[p].append(x)
        h = h.group("process", grp_map)

        # Remove the categories which were already correlated together so as to not double count
        if already_correlated:
            for k in already_correlated:
                if self.verbose: print(f"Removing: {k}")
            h = h.remove("systematic", list(already_correlated))

        return h

    def correlate_years_regular(self, h):
        """Sums up years for a regular hist.Hist object (like LHEPDFweights)."""
        # Identify common process names
        all_procs = list(h.axes["process"])
        unique_procs = set([re.sub(r'UL\d+(APV)?$', '', p) for p in all_procs])
        
        new_proc_axis = hist.axis.StrCategory(unique_procs, name="process", growth=True)
        other_axes = [ax for ax in h.axes if ax.name != "process"]
        h_corr = hist.Hist(new_proc_axis, *other_axes, storage=h.storage_type())
        
        for p_base in unique_procs:
            # Match all years for this process
            matching_olds = [p for p in all_procs if re.match(rf"^{p_base}_UL\d+(APV)?$", p)]
            if not matching_olds: continue
            
            summed = h[{"process": matching_olds}][{"process": sum}]
            h_corr[{"process": p_base}] = summed.view(flow=True)
            
        return h_corr
    
    # def add_pdf_shapes(self, km_dist):
    #     """Calculates PDF scaling by re-mapping the histogram to include new categories."""
    #     if "LHEPDFweights" not in self.hists:
    #         return

    #     h_main = self.hists[km_dist]
    #     h_pdf = self.hists["LHEPDFweights"]

    #     # 1. Create a new axis that includes the original categories plus the PDF ones
    #     old_syst_axis = h_main.axes["systematic"]
    #     new_syst_list = list(old_syst_axis) + ["PDFUp", "PDFDown"]
    #     new_syst_axis = hist.axis.StrCategory(new_syst_list, name="systematic", growth=True)

    #     # 2. Create a new histogram with the expanded axis
    #     # We keep all other axes the same as h_main
    #     other_axes = [ax for ax in h_main.axes if ax.name != "systematic"]
    #     h_new = hist.Hist(new_syst_axis, *other_axes, storage=h_main.storage_type())

    #     # 3. Copy existing data from h_main to h_new
    #     for syst in old_syst_axis:
    #         h_new[{"systematic": syst}] = h_main[{"systematic": syst}].view(flow=True)

    #     # Replace the old histogram in the dictionary
    #     self.hists[km_dist] = h_new
    #     h_main = self.hists[km_dist] # Point h_main to the new expanded object

    #     # Now define the indices for the loop
    #     idx_nom = h_main.axes["systematic"].index("nominal")
    #     idx_up  = h_main.axes["systematic"].index("PDFUp")
    #     idx_dn  = h_main.axes["systematic"].index("PDFDown")

    #     view = h_main.view(flow=True)

    #     for ch in h_main.axes["channel"]:
    #         ch_idx = h_main.axes["channel"].index(ch)
    #         for p in h_main.axes["process"]:
    #             p_idx = h_main.axes["process"].index(p)
    #             if p not in h_pdf.axes["process"]:
    #                 continue

    #             # PDF calculation logic (Standard Hessian/MC replicas)
    #             pdf_vals = h_pdf[{"channel": ch, "process": p}].values()
    #             nom_pdf = pdf_vals[..., 0]
    #             replicas = pdf_vals[..., 1:]
    #             sigma_pdf = np.sqrt(np.sum(np.square(replicas - nom_pdf[..., np.newaxis]), axis=-1))

    #             # Scaling Factors
    #             rel_sigma = np.divide(sigma_pdf, nom_pdf, out=np.zeros_like(sigma_pdf), where=nom_pdf!=0)
    #             scale_up = 1.0 + rel_sigma
    #             scale_dn = np.maximum(1.0 - rel_sigma, 0.0)

    #             # Apply to all EFT coefficients in the HistEFT storage
    #             for wc_key, coeff_array in view.items():
    #                 nom_vals = coeff_array[ch_idx, p_idx, idx_nom]['value']
    #                 nom_vars = coeff_array[ch_idx, p_idx, idx_nom]['variance']

    #                 # Fill the new PDF slots (1:-1 excludes under/overflow)
    #                 coeff_array[ch_idx, p_idx, idx_up, 1:-1]['value'] = nom_vals[1:-1] * scale_up
    #                 coeff_array[ch_idx, p_idx, idx_up, 1:-1]['variance'] = nom_vars[1:-1]
                    
    #                 coeff_array[ch_idx, p_idx, idx_dn, 1:-1]['value'] = nom_vals[1:-1] * scale_dn
    #                 coeff_array[ch_idx, p_idx, idx_dn, 1:-1]['variance'] = nom_vars[1:-1]


    # def add_pdf_shapes(self, km_dist):
    #     """Calculates PDF and injects into growth-enabled systematic axis."""
    #     h_main = self.hists[km_dist]
    #     h_pdf = self.hists["LHEPDFweights"]

    #     # Note: h_main is grouped/year-correlated, h_pdf is now grouped/year-correlated
    #     for ch in h_main.axes["channel"]:
    #         for p in h_main.axes["process"]:
    #             if p not in h_pdf.axes["process"]:
    #                 print(f"DEBUG: Process {p} not found in PDF histogram!")
    #                 continue
    #             print(f"DEBUG: Adding PDF shapes for {p}")

    #             # Get PDF info from LHEPDFweights
    #             # PDFindex 0 = nominal, 1-103 = replicas
    #             pdf_vals = h_pdf[{"channel": ch, "process": p}].values()
    #             nom_pdf = pdf_vals[..., 0]
    #             PDFvariations = pdf_vals[..., 1:]

    #             # Calculate standard deviation (Hessian/MC replicas)
    #             diff_sq = np.square(PDFvariations - nom_pdf[..., np.newaxis])
    #             sigma_pdf = np.sqrt(np.sum(diff_sq, axis=-1))

    #             # Get the nominal yield from the kinematic distribution (mllbb)
    #             h_nom_main = h_main[{"channel": ch, "process": p, "systematic": "nominal"}]
    #             nom_main_vals = h_nom_main.values()
    #             nom_main_vars = h_nom_main.variances()

    #             # Calculate Up/Down: nom +/- sigma
    #             # Because it's correlated, we apply the absolute sigma 
    #             # calculated from the PDF-dedicated histogram
    #             up_vals = nom_main_vals + sigma_pdf
    #             down_vals = np.maximum(nom_main_vals - sigma_pdf, 0)

    #             # Inject into main histogram (Growth axis handles the new strings)
    #             h_main.view(flow=True)[{"channel": ch, "process": p, "systematic": "PDFUp"}] = np.stack([up_vals, nom_main_vars], axis=-1)
    #             h_main.view(flow=True)[{"channel": ch, "process": p, "systematic": "PDFDown"}] = np.stack([down_vals, nom_main_vars], axis=-1)


    def channels(self, km_dist):
        return list(self.hists[km_dist].axes["channel"])

    def processes(self, km_dist):
        return list(self.hists[km_dist].axes["process"])
    
    def make_scalings_json(self,scalings_json,ch,km_dist,p,wc_names,scalings):
        print(f"making scalings json for {km_dist}, {p}")
        scalings = scalings.tolist()
        scalings_json.append(
            {
                "channel": ch + "_" + str(km_dist),
                "process": p + "_sm",  # NOTE: needs to be in the datacard
                "parameters": ["cSM[1]"] + [
                    self.format_wc(wcname) for wcname in wc_names
                ],
                "scaling":
                    scalings[1:], # exclude underflow bin
            }
        )
        return scalings_json

    def format_wc(self,wcname):
        if wcname in self.rotate:
            return self.rotate[wcname]
        lo, hi = self.wc_ranges[wcname]
        return "%s[0,%.1f,%.1f]" % (wcname, lo, hi)
    
    def get_selected_wcs(self,km_dist,ch_lst=[]):
        """
            For each process, iterates over every channel and every bin checking the EFT parameterization
            coefficients for if they have a significant impact or not relative to the SM contribution. If
            any term from any channel+bin is determined to be significant, the WC is selected, otherwise
            it is excluded for that process and won't be included in the EFT decomposition
        """
        tic = time.time()
        h = self.hists[km_dist].integrate("systematic",["nominal"])
        if ch_lst:
            # Only select from a subset of channels
            if self.verbose:
                print(f"Selecting WCs from subset of channels: {ch_lst}")
            h.prune("channel", ch_lst)

        procs = list(h.axes["process"])
        selected_wcs = {p: set() for p in procs}

        wcs = ["sm"] + h.wc_names

        # This maps a WC to a list whose elements are the indices of the coefficient array of the
        #   HistEFT that involve that particular WC
        # NOTE: Building up the index array MUST match exactly with how the HistEFT coeff array is
        #       constructed/computed [1], otherwise the index array that gets computed won't pick
        #       out the correct coeff array indices for the corresponding WC!
        # [1] https://github.com/TopEFT/topcoffea/blob/3bef686fead216183ebb6dfb464e67629cfe75f5/topcoffea/modules/eft_helper.py#L32-L36
        wc_to_terms = {}
        start_index = 1
        index = start_index
        for i in range(len(wcs)):
            wc1 = wcs[i]
            wc_to_terms[wc1] = set()
            for j in range(i+1):
                wc2 = wcs[j]
                wc_to_terms[wc1].add(index)
                wc_to_terms[wc2].add(index)
                index += 1

        # Convert the set to a sorted np.array
        for wc in wcs:
            wc_to_terms[wc] = np.array(sorted(wc_to_terms[wc]))

        for p in procs:
            if not self.is_signal(p):
                continue
            p_hist = h.integrate("process",[p])
            for wc,idx_arr in wc_to_terms.items():
                if len(self.coeffs) and not wc in self.coeffs:
                    continue
                if wc == "sm":
                    continue
                if wc == "ctlTi" and p == "tttt":
                    continue
                for sp_key, arr in p_hist.view(as_dict=True, flow=True).items():
                    # Ignore underflow, and overflow bins
                    sl_arr = arr[1:-1]
                    # Here we replace any SM terms that are too small with a large dummy value
                    sm_norm = np.where(sl_arr[:,start_index] < 1e-5,999,sl_arr[:,start_index])
                    # Normalize each sub-array by corresponding SM term
                    n_arr = (sl_arr.T / sm_norm).T
                    # This will contain only the coefficients which involve the given WC
                    wc_terms = np.abs(n_arr[:,idx_arr])
                    if np.any(wc_terms > self.tolerance):
                        selected_wcs[p].add(wc)
                        break
        if self.verbose:
            dt = time.time() - tic
            print(f"WC Selection Time: {dt:.2f} s")
        return selected_wcs
    
    def analyze(self,km_dist,ch,selected_wcs, crop_negative_bins, wcs_dict):
        """ Handles the EFT decomposition and the actual writing of the ROOT and text datacard files."""
        if not km_dist in self.hists:
            print(f"[ERROR] Unknown kinematic distribution: {km_dist}")
            return None
        elif ch not in self.hists[km_dist].axes["channel"]:
            print(f"[ERROR] Unknown channel {ch}")
            return None

        print(f"Analyzing {km_dist} in {ch}")

        bin_str = f"bin_{ch}_{km_dist}"
        col_width = max(PRECISION*2+5,len(bin_str))
        syst_width = 0

        if km_dist != "njets":
            num_j,num_b = self.get_jet_mults(ch)
        else:
            num_j,num_b = 0,0
        num_l = self.get_lep_mult(ch)

        outf_root_name = self.FNAME_TEMPLATE.format(cat=ch,kmvar=km_dist,ext="root")

        h = self.hists[km_dist]
        h_sumw2 = None

        if f"{km_dist}_sumw2" in self.hists:
            h_sumw2 = self.hists[km_dist+"_sumw2"]
        elif ("systematic" in h.axes.name) and ("sumw2" in h.axes["systematic"]):
            h_sumw2 = h[{"systematic": "sumw2"}]
            print(f"sumw2 found in systematic axis")
        else:
            print("No sumw2 histogram found! Setting errors to 0")

        ch_hist = h.integrate("channel",[ch])
        print(f"ch_hist:{ch_hist}")
        ch_sumw2 = h_sumw2.integrate("channel", [ch]) if h_sumw2 is not None else None
        data_obs = np.zeros((2, ch_hist.dense_axis.extent))
            
        print(f"Generating root file: {outf_root_name}")
        tic = time.time()
        num_h = 0
        all_shapes = set()
        text_card_info = {}
        outf_root_name = os.path.join(self.out_dir,outf_root_name)
        with uproot.recreate(outf_root_name) as f:
            for p,wcs in selected_wcs.items():

                # --- START OF PDF BLOCK ---
                # Check if we should calculate PDF shapes for this process
                is_pdf_skipped = p in self.syst_skip.get("PDF", [])
                if not is_pdf_skipped and "LHEPDFweights" in self.hists:
                # if "PDF" not in self.syst_skip.get("PDF", []) or p not in self.syst_skip.get("PDF", []):
                    # if "LHEPDFweights" in self.hists and p in self.hists["LHEPDFweights"].axes["process"]:
                    if p in self.hists["LHEPDFweights"].axes["process"]:

                        ### this version worked but didn't include overflow ### 
                        # nom_hist_vals = (ch_hist[{'process':p, 'systematic':'nominal'}].as_hist({}).values(flow=True))[0]
                        # nom_hist_vals = (ch_hist[{'process':p, 'systematic':'nominal'}].as_hist({}).values())

                        # h_pdf = self.hists["LHEPDFweights"]
                        # h_pdf = h_pdf[{"channel": ch, "process": p}]

                        # # print(f"nominal hist: { ch_hist[{'process':p, 'systematic':'nominal'}]}")
                        # print(f"nominal hist vals: {nom_hist_vals}")
                        # # print(f"h_pdf: {h_pdf}")

                        # nominal = h_pdf[{'PDFindex':0}].values() #flow=True
                        # variations = h_pdf[{"PDFindex": slice(1, None)}].values() #flow=True

                        # diff_sq = np.square(variations - nominal[..., np.newaxis])
                        # sigma_pdf = np.sqrt(np.sum(diff_sq, axis=-1)) # [1:]

                        # print(f"sigma_PDF: {sigma_pdf}")

                        # pdf_up = nom_hist_vals + sigma_pdf
                        # pdf_down = np.maximum(nom_hist_vals - sigma_pdf, 0)

                        # # pdf_var = h_pdf[{"channel": ch, "process": p}].variances()[..., 0]
                        # pdf_var = h_pdf.variances()[..., 0]

                        # # np.append(pdf_up, (nom_hist_vals[-1]*1.015))
                        # # np.append(pdf_down, (nom_hist_vals[-1] - (nom_hist_vals[-1]*0.015)))
                        # print(f"pdf_up: {pdf_up}")
                        # print(f"pdf_down: {pdf_down}")
                        # print(f"pdf_var: {pdf_var}")

                        ### trying to get the overflow bin
                        nom_hist_vals = ch_hist[{'process':p, 'systematic':'nominal'}].as_hist({}).values()[0]

                        h_pdf = self.hists["LHEPDFweights"]
                        h_pdf = h_pdf[{"channel": ch, "process": p}]

                        nominal = h_pdf[{'PDFindex':0}].values() #flow=True
                        variations = h_pdf[{"PDFindex": slice(1, None)}].values() #flow=True
                        diff_sq = np.square(variations - nominal[..., np.newaxis])
                        sigma_pdf = np.sqrt(np.sum(diff_sq, axis=-1))

                        nom_overflow = ch_hist[{'process':p, 'systematic':'nominal'}].as_hist({}).values(flow=True)[0][-1]
                        pdf_nom_overflow = h_pdf[{'PDFindex': 0}].values(flow=True)[-1]
                        pdf_vars_overflow = h_pdf[{"PDFindex": slice(1, None)}].values(flow=True)[:, -1]
                        sigma_overflow = np.sqrt(np.sum(np.square(pdf_vars_overflow - pdf_nom_overflow)))

                        nom_full = np.append(nom_hist_vals, nom_overflow)
                        sigma_full = np.append(sigma_pdf, sigma_overflow)
                        # pdf_var_full = h_pdf.variances(flow=True)[..., 0][1:]
                        pdf_var_full = h_pdf[{"PDFindex": 0}].variances(flow=True)[1:]

                        rel_sigma = sigma_full / nom_full
                        # If the overflow relative uncertainty is > 50%, cap it at the last interior bin's uncertainty
                        if rel_sigma[-1] > 0.5: 
                            print(f"Warning: Overflow PDF uncertainty ({rel_sigma[-1]:.2%}) is suspicious. Capping.")
                            sigma_full[-1] = nom_full[-1] * rel_sigma[-2]

                        print(f"nom_full: {nom_full}")
                        print(f"sigma_PDF: {sigma_full}")
                        pdf_up = nom_full + sigma_full
                        pdf_down = np.maximum(nom_full - sigma_full, 0)
                        print(f"pdf_up: {pdf_up}")
                        print(f"pdf_down: {pdf_down}")
                        print(f"pdf_var: {pdf_var_full}")

                        for syst_name, vals in [("PDFUp", pdf_up), ("PDFDown", pdf_down)]:
                            hist_name = f"{p}_sm_{syst_name}"
                            arr = [vals, pdf_var_full] 
                            # f[hist_name] = to_hist(arr, hist_name, zero_wgts=(p != "fakes"))
                            f[hist_name] = to_hist(arr, hist_name, zero_wgts=None)
                            
                            # Ensure the process entry exists in the dictionary
                            if f"{p}_sm" not in text_card_info:
                                text_card_info[f"{p}_sm"] = {"shapes": set(), "rate": -1}
                            
                            text_card_info[f"{p}_sm"]["shapes"].add("PDF")
                            all_shapes.add("PDF")

                # --- END OF PDF BLOCK ---

                proc_hist = ch_hist.integrate("process",[p])
                if ch_sumw2 is not None and p in ch_sumw2.axes["process"]:
                    proc_sumw2 = ch_sumw2.integrate("process", [p])
                else:
                    # If the process isn't in sumw2, it has no statistical variance entries
                    proc_sumw2 = None
                    if self.verbose:
                        print(f"Note: Process {p} has no sumw2 entries (likely zero yield).")

                # proc_sumw2 = ch_sumw2 if ch_sumw2 is None else ch_sumw2.integrate("process",[p])
                if self.verbose:
                    print(f"Decomposing {ch}-{p}")
                decomposed_templates = self.decompose(proc_hist,proc_sumw2,wcs)
                is_eft = self.is_signal(p)
                # Note: This feels like a messy way of picking out the data_obs info
                if p == "data":
                    data_sm = decomposed_templates.pop("sm")
                    if self.use_real_data:
                        print(f"getting real data!")
                        if len(data_sm) != 1:
                            data_sm = next((arr[0] for key, arr in data_sm.items() if key.systematic == 'nominal'), None)
                            print(data_sm)
                        if data_sm is None:
                            raise RuntimeError("Nominal data not found in templates")
                        if len(data_sm) == 0:
                            raise RuntimeError("obs data has unexpected number of sparse bins")
                            # if len(data_sm) != 1:
                            #     raise RuntimeError("obs data has unexpected number of sparse bins")
                        elif sum(data_obs[0]) != 0:
                            raise RuntimeError("filling obs data more than once!")
                        data_obs += data_sm
                        # for sp_key,arr in data_sm.items():
                        #     data_obs += arr
                if not self.use_AAC:
                    decomposed_templates = {k: v for k, v in decomposed_templates.items() if k == 'sm'}
                for base,v in decomposed_templates.items():
                    proc_name = f"{p}_{base}"
                    col_width = max(len(proc_name),col_width)
                    if proc_name not in text_card_info:
                        text_card_info[proc_name] = {
                            "shapes": set(),
                            "rate": -1
                        }
                    # There should be only 1 sparse axis at this point, the systematics axis
                    check_zero_arr0 = False
                    check_zero_arr1 = False

                    for sp_key,arr in v.items():
                        syst = sp_key.systematic

                        if syst == "sumw2": continue

                        syst_base = self.strip_fluctuation(syst)
                        if syst_base in self.syst_skip:
                            if p in self.syst_skip[syst_base]:
                                if self.verbose:
                                    print(f"\tSkipping shape systematic {syst_base} for process {p}")
                                continue

                        if crop_negative_bins:
                            negative_bin_mask = np.where( arr[0] < 0) # see where bins are negative
                            arr[0][negative_bin_mask] = np.zeros_like( arr[0][negative_bin_mask] )  # set those to zero
                            if arr[1] is not None:
                                arr[1][negative_bin_mask] = np.zeros_like( arr[1][negative_bin_mask] )  # if there's a sumw2 defined, that one's set to zero as well. Otherwise we will get 0 +/- something, which is compatible with negative

                        if syst =="nominal":  # check systematics error for fake factors
                            if sum(arr[0]) == 0:
                                check_zero_arr0 = True
                            if sum(arr[1]) == 0:
                                check_zero_arr1 = True
                        if "FF" in syst:
                            if check_zero_arr0 and sum(arr[0]) != 0:
                                raise Warning("Systematics Error arr[0]:Zero values in 'nominal' but non-zero in '%s'" % (syst))
                            if check_zero_arr1 and sum(arr[1]) != 0:
                                raise Warning("Systematics Error arr[1]:Zero values in 'nominal' but non-zero in '%s'" % (syst))

                        sum_arr = sum(arr[0])
                        if sum_arr == 0: continue #TODO find a more elegant solution
                        if syst == "nominal" and base == "sm":
                            if self.verbose:
                                print(f"\t{proc_name:<12}: {sum_arr:.4f} {arr[0]}")
                            if not self.use_real_data:
                                # Create asimov dataset
                                vals = wcs_dict # set wcs to certain values from command line
                                decomposed_templates_Asimov = self.decompose(proc_hist,proc_sumw2,wcs,vals)
                                data_sm = decomposed_templates_Asimov.pop("sm")
                                data_obs += data_sm[sp_key]
                        if syst == "nominal":
                            hist_name = f"{proc_name}"
                            text_card_info[proc_name]["rate"] = sum_arr
                        else:
                            hist_name = f"{proc_name}_{syst}"
                            
                            # Systematics in the text datacard don't have the Up/Down postfix
                            syst_base = self.strip_fluctuation(syst) #replace("Up","").replace("Down","")
                            if syst_base in self.syst_shape_decorrelate:
                                # We want to split this systematic to be uncorrelated between certain
                                #   processes, so we modify the systematic name to make combine treat
                                #   them as separate systematics. Also, we use 'p' instead of 'proc_name'
                                #   for renaming since we want the decomposed EFT terms for a particular
                                #   process to share the same nuisance parameter
                                matched = []
                                for r in self.syst_shape_decorrelate[syst_base]:
                                    if regex_match([p],r["matches"]):
                                        # The matched process should have this systematic put into a new group
                                        matched.append(r["group"])
                                if len(matched) == 0:
                                    # No matches found, so keep the original systematic name
                                    split_syst = syst_base
                                elif len(matched) == 1:
                                    # Found a match, so decorrelate the process from non-matched processes
                                    group = matched[0]
                                    split_syst = f"{syst_base}_{group}"
                                    if group == "":
                                        # In the special case that the group is an empty string,
                                        #   decorrelate ALL matched processes
                                        split_syst = f"{syst_base}_{p}"
                                else:
                                    # We shouldn't have more than one match for a given systematic
                                    raise RuntimeError(f"Unable to decorrelate shape systematic {syst_base} for {p}. Multiple group matches found: {matched}")
                                hist_name = hist_name.replace(syst_base,split_syst)
                                all_shapes.add(split_syst)
                                text_card_info[proc_name]["shapes"].add(split_syst)
                                if base == "sm" and self.verbose:
                                    print(f"\tDecorrelate {p} for {syst_base} into {split_syst} ({syst.replace(syst_base,'')})")
                            else:
                                all_shapes.add(syst_base)
                                text_card_info[proc_name]["shapes"].add(syst_base)
                            syst_width = max(len(syst),syst_width)
                        zero_out_sumw2 = False # p != "fakes" and "close" not in p # Zero out sumw2 for all proc but fakes, so that we only do auto stats for fakes
                        f[hist_name] = to_hist(arr,hist_name,zero_wgts=zero_out_sumw2)

                        num_h += 1

                # obtain the scalings for scalings.json file
                print(f"p: {p}")
                print(f"self.SIGNALS: {self.SIGNALS}")
                if p in self.SIGNALS:
                    if self.wc_scalings:
                        # scalings = h[{'channel':ch,'process':p,'systematic':'nominal'}].make_scaling(flow='show', wc_list=self.wc_scalings)
                        scalings = h[{'channel':ch,'process':p,'systematic':'nominal'}].make_scaling(flow='show')
                        self.scalings_json = self.make_scalings_json(self.scalings,ch,km_dist,p,self.wc_scalings,scalings)
                    else:
                        scalings = h[{'channel':ch,'process':p,'systematic':'nominal'}].make_scaling(flow='show')
                        self.scalings_json = self.make_scalings_json(self.scalings,ch,km_dist,p,h.wc_names,scalings)
            f["data_obs"] = to_hist(data_obs,"data_obs")

        line_break = "##----------------------------------\n"
        left_width = len(line_break) + 2
        left_width = max(syst_width+len("shape")+1,left_width)

        outf_card_name = self.FNAME_TEMPLATE.format(cat=ch,kmvar=km_dist,ext="txt")
        print(f"Generating text file: {outf_card_name}")
        outf_card_name = os.path.join(self.out_dir,outf_card_name)
        with open(outf_card_name,"w") as f:
            f.write(f"shapes *        * {os.path.split(outf_root_name)[1]} $PROCESS $PROCESS_$SYSTEMATIC\n")
            f.write(line_break)
            f.write(f"bin         {bin_str}\n")
            f.write(f"observation {sum(data_obs[0]):.{PRECISION}f}\n")
            f.write(line_break)
            f.write(line_break)

            # Note: This list is what controls the columns in the text datacard, if a process appears
            #       in this list it should NEVER be skipped in any of the following for loops.
            # proc_order = sorted(text_card_info.keys())
            proc_order = [k for k in text_card_info.keys() if text_card_info[k]["rate"] != -1]  # rate = -1 only happens when there's no syst histograms (e.g. flips in 3l/4l)

            # Bin row
            row = [f"{'bin':<{left_width}}"]    # Python string formatting is pretty great!
            for p in proc_order:
                row.append(f"{bin_str:>{col_width}}")
            row = " ".join(row) + "\n"
            f.write(row)

            # 1st process row
            row = [f"{'process':<{left_width}}"]
            for p in proc_order:
                row.append(f"{p:>{col_width}}")
            row = " ".join(row) + "\n"
            f.write(row)

            # 2nd process row
            row = [f"{'process':<{left_width}}"]
            bkgd_count =  1
            sgnl_count = -1
            for p in proc_order:
                if any([x in p for x in self.SIGNALS]): # Check for if the process is signal or not
                    row.append(f"{sgnl_count:>{col_width}}")
                    sgnl_count += -1
                else:
                    row.append(f"{bkgd_count:>{col_width}}")
                    bkgd_count += 1
            row = " ".join(row) + "\n"
            f.write(row)

            # Rate row
            row = [f"{'rate':<{left_width}}"]
            for p in proc_order:
                r = text_card_info[p]["rate"]
                if r < 0:
                    print(f"Process {p} has negative total rate: {r:.3f} -> setting to 0 in text card")
                    r = 0
                row.append(f"{r:>{col_width}.{PRECISION}f}") # Do not challenge me on Python string formatting!
            row = " ".join(row) + "\n"
            f.write(row)
            f.write(line_break)

            # Shape systematics rows
            for syst in sorted(all_shapes):
                left_text = f"{syst:<{syst_width}} shape"
                row = [f"{left_text:<{left_width}}"]
                for p in proc_order:
                    if syst in text_card_info[p]["shapes"]:
                        row.append(f"{'1':>{col_width}}")
                    else:
                        row.append(f"{'-':>{col_width}}")
                row = " ".join(row) + "\n"
                f.write(row)

            # Rate systematics rows
            for k,rate_syst in self.rate_systs.items():
                syst_name = rate_syst.name
                left_text = f"{syst_name:<{syst_width}} lnN"
                row = [f"{left_text:<{left_width}}"]
                for p in proc_order:
                    proc_name = self.get_process(p) # Strips off any "_sm" or "_lin_*" junk
                    v = rate_syst.get_process(proc_name)
                    row.append(f"{v:>{col_width}}")
                row = " ".join(row) + "\n"
                f.write(row)

            if self.do_mc_stat:
                f.write("* autoMCStats 10\n")
            else:
                f.write("* autoMCStats -1\n")

        outf_json_name = self.FNAME_TEMPLATE.format(cat=ch,kmvar=km_dist,ext="json")
        with open(os.path.join(self.out_dir,f"{outf_json_name}"),"w") as f:
            print('making scalings json:', os.path.join(self.out_dir,f"{outf_json_name}"))
            json.dump(self.scalings_json, f, indent=4)

        dt = time.time() - tic
        print(f"File Write Time: {dt:.2f} s")
        print(f"Total Hists Written: {num_h}")
                        
    def decompose(self,h,sumw2,wcs,vals={}):
        """
            Decomposes the EFT quadratic parameterization coefficients into combinations that result
            in non-negative coefficient terms.

            Note: All other WCs are assumed set to 0
            sm piece:    set(c1=0.0)
            lin piece:   set(c1=1.0)
            mixed piece: set(c1=1.0,c2=1.0)
            quad piece:  0.5*[set(c1=2.0) - 2*set(c1=1.0) + set(sm)]
        """
        tic = time.time()

        sm = h.eval({})

        if vals:
            # If we are evaluating at an EFT point, but sumw2 has no EFT info
            try:
                sm_w2 = sumw2.eval(vals) if sumw2 is not None else None
            except (LookupError, KeyError): #
                sm_w2 = sumw2.eval({}) 
            
        else:
            # sm_w2 = sumw2.eval({})
            sm_w2 = sumw2.eval({}) if sumw2 is not None else None

        sm = add_sumw2_stub(sm,sm_w2)

        # Note: The keys of this dictionary are a pretty contrived, but are useful later on
        r = {}
        r["sm"] = sm
        terms = 1
        for n1, wc1 in enumerate(wcs):
            tmp_lin_1 = h.eval({wc1: 1.0})
            tmp_lin_2 = h.eval({wc1: 2.0})

            tmp_lin_1 = add_sumw2_stub(tmp_lin_1)
            tmp_lin_2 = add_sumw2_stub(tmp_lin_2)

            lin_name = f"lin_{wc1}"
            quad_name = f"quad_{wc1}"

            terms += 2

            r[lin_name] = tmp_lin_1
            r[quad_name] = {}
            for sp_key in h.categorical_keys:
                r[quad_name][sp_key] = []
                for i in range(2):
                    r[quad_name][sp_key].append(
                        0.5 * (
                            tmp_lin_2[sp_key][i] - 2 * tmp_lin_1[sp_key][i] + sm[sp_key][i]
                        )
                    )

            for n2, wc2 in enumerate(wcs):
                if n1 >= n2:
                    continue
                mixed_name = f"quad_mixed_{wc1}_{wc2}"
                mixed = h.eval({wc1: 1.0, wc2: 1.0})
                mixed = add_sumw2_stub(mixed)
                r[mixed_name] = mixed
                terms += 1

        toc = time.time()
        dt = toc - tic
        if self.verbose:
            print(f"\tDecompose Time: {dt:.2f} s")
            print(f"\tTotal terms: {terms}")

        return r
