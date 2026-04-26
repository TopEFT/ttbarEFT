import os
import time
import json
import shutil
import argparse
import numpy as np
import functools, operator

from topcoffea.modules.utils import regex_match,clean_dir,dict_comp
from ttbarEFT.modules.new_datacard_tools import *


def combine_vineReduce_channels(pkl): 
    vineReduce_dict = load_pickle(pkl)
    first_ch = next(iter(vineReduce_dict))
    variables = vineReduce_dict[first_ch].keys()
    
    hists = {}
    for var in variables: 
        h_list = [vineReduce_dict[ch][var] for ch in vineReduce_dict]
        hists[var] = functools.reduce(operator.add, h_list)

    return hists

def flatten_vineReduce_hists(pkl, outname = 'flat_hists'):
    # Flattens the histogram structure from vineReduce and saves as a new pickle. 
    # Returns the name of the new pickle file 
    
    flattend_hists = combine_vineReduce_channels(pkl)
    tmp_pkl = f"{outname}.pkl.gz"
    with gzip.open(tmp_pkl, "wb") as f:
        pickle.dump(flattend_hists, f)
    
    return tmp_pkl

def load_pickle(fname):
    return pickle.load(gzip.open(fname))


def run_local(dc,km_dists,channels,selected_wcs, crop_negative_bins, wcs_dict):
    for km_dist in km_dists:
        all_chs = dc.channels(km_dist)
        matched_chs = regex_match(all_chs,channels)
        if channels:
            print(f"Channels to process: {matched_chs}")
        for ch in matched_chs:
            r = dc.analyze(km_dist,ch,selected_wcs, crop_negative_bins, wcs_dict)


def main():
    parser = argparse.ArgumentParser(description="You can select which file to run over")
    parser.add_argument("pkl_file",nargs="?",help="Pickle file with histograms to run over")
    parser.add_argument("--rate-syst-json","-s",default="params/rate_systs.json",help="Rate related systematics json file, path relative to topeft_path()")
    parser.add_argument("--miss-parton-file","-m",default="data/missing_parton/missing_parton.root",help="File for missing parton systematic, path relative to topeft_path()")
    parser.add_argument("--selected-wcs-ref",default="test/selectedWCs.json",help="Reference file for selected wcs")
    parser.add_argument("--out-dir","-d",default=".",help="Output directory to write root and text datacard files to")
    parser.add_argument("--var-lst",default=['mllbb'],action="extend",nargs="+",help="Specify a list of variables to make cards for.")
    parser.add_argument("--ch-lst","-c",default=[],action="extend",nargs="+",help="Specify a list of channels to process.")
    parser.add_argument("--do-mc-stat",action="store_true",help="Add bin-by-bin statistical uncertainties with the autoMCstats option (for background)")
    parser.add_argument("--ignore","-i",default=[],action="extend",nargs="+",help="Specify a list of processes to exclude, must match name from 'sample' axis modulo UL year")
    parser.add_argument("--drop-syst",default=[],action="extend",nargs="+",help="Specify one or more template systematics to remove from the datacard")
    parser.add_argument("--POI",default=[],help="List of WCs (comma separated)")
    parser.add_argument("--year","-y",default=[],action="extend",nargs="+",help="Run over a subset of years")
    parser.add_argument("--do-nuisance",action="store_true", help="Include nuisance parameters")
    parser.add_argument("--unblind",action="store_true",help="If set, use real data, otherwise use asimov data")
    parser.add_argument("--verbose","-v",action="store_true",help="Set to verbose output")
    parser.add_argument("--select-only",action="store_true",help="Only run the WC selection step")
    parser.add_argument("--skip-selected-wcs-check",action="store_true",help="Do not raise an error if the selected WCs disagree with ref")
    parser.add_argument("--use-selected",default="",help="Load selected process+WC combs from a file. Skips doing the normal selection step.")
    parser.add_argument("--condor","-C",action="store_true",help="Split up the channels into multiple condor jobs")
    parser.add_argument("--chunks","-n",default=1,help="The number of channels each condor job should process")
    parser.add_argument("--keep-negative-bins",action="store_true",help="Don't crop negative bins")
    parser.add_argument("--use-AAC","-A",action="store_true",help="Include all EFT templates in datacards for AAC model")
    parser.add_argument("--wc-vals", default="",action="store", nargs="+", help="Specify the corresponding wc values to set for the wc list")
    parser.add_argument("--wc-scalings", default=[],action="extend",nargs="+",help="Specify a list of wc ordering for scalings.json")

    args = parser.parse_args()
    pkl_file   = args.pkl_file
    rs_json    = args.rate_syst_json
    mp_file    = args.miss_parton_file
    out_dir    = args.out_dir
    years      = args.year
    var_lst    = args.var_lst
    ch_lst     = args.ch_lst
    do_mc_stat = args.do_mc_stat
    wcs        = args.POI
    ignore     = args.ignore
    do_nuis    = args.do_nuisance
    drop_syst  = args.drop_syst
    unblind    = args.unblind
    verbose    = args.verbose
    use_AAC     = args.use_AAC
    wc_vals    = args.wc_vals

    wc_scalings = args.wc_scalings
    select_only = args.select_only
    use_selected = args.use_selected

    use_condor = args.condor
    chunks = int(args.chunks)

    if isinstance(wcs,str):
        wcs = wcs.split(",")

    kwargs = {
        "wcs": wcs,
        "rate_syst_path": rs_json,
        # "missing_parton_path": mp_file,
        "out_dir": out_dir,
        "var_lst": var_lst,
        "do_mc_stat": do_mc_stat,
        "ignore": ignore,
        "do_nuisance": do_nuis,
        "drop_syst": drop_syst,
        "unblind": unblind,
        "verbose": verbose,
        "year_lst": years,
        "use_AAC":  use_AAC,
        "wc_vals": wc_vals,
        "wc_scalings": wc_scalings,
    }

    prepared_pkl = flatten_vineReduce_hists(pkl_file)

    if out_dir != "." and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Copy over make_cards.py ASAP so a user can't accidentally modify it before the submit jobs run
    if use_condor and not os.path.samefile(os.getcwd(),out_dir):
        shutil.copy("make_cards.py",out_dir)

    tic = time.time()
    # dc = DatacardMaker(pkl_file,**kwargs)
    dc = DatacardMaker(
        prepared_pkl,
        year_lst=[],         # Processes all years found in pkl
        do_nuisance=True,    # Set to False if you want to skip systematics entirely for now
        do_mc_stat=True,
        out_dir=out_dir
    )

    dists = var_lst if len(var_lst) else dc.hists.keys()

    ch_lst = []
    selected_wcs = {}
    for km_dist in dists:
        all_chs = dc.channels(km_dist)
        matched_chs = regex_match(all_chs,ch_lst)
        
        dist_wcs = dc.get_selected_wcs(km_dist,matched_chs)
        for p,wcs in dist_wcs.items():
            if not p in selected_wcs:
                selected_wcs[p] = []
            for wc in wcs:
                if not wc in selected_wcs[p]:
                    selected_wcs[p].append(wc)
                    
    with open(os.path.join(out_dir,"selectedWCs.txt"),"w") as f:
        selected_wcs_for_json = {}
        for p,v in selected_wcs.items():
            if not dc.is_signal(p):
                # WC selection will include backgrounds in the dict (always with an empty list), so remove them here
                continue
            selected_wcs_for_json[p] = list(v)
        json.dump(selected_wcs_for_json,f)

    # set wcs_dict to be zero - we want the asimov dataset made with all WC set to zero
    wc_vals = ""
    wc_vals = ''.join(wc_vals)
    wcs_dict = eval("dict({})".format(wc_vals))


    run_local(dc,dists,ch_lst,selected_wcs, not args.keep_negative_bins, wcs_dict)

    # make pre-selection scalings.json
    print("Making scalings-preselect.json file...")
    with open(os.path.join(out_dir,"scalings-preselect.json"),"w") as f:
        json.dump(dc.scalings, f, indent=4)

    dt = time.time() - tic
    print(f"Total Time: {dt:.2f} s")
    print("Finished!")

main()
