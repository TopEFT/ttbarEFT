import os
import json
import argparse
import topcoffea.modules.utils as utils

def make_sample_json(xsec, year, treeName, histAxisName, options, era, paths, outname):

    # check that output json file doesn't already exist
    outputFile = outname+'.json'
    if os.path.exists(outputFile):
        assert not os.path.exists(outputFile), "Output file will not be overwritten! Backup the file and try again."

    # save input arguments
    sampdic = {}
    sampdic['xsec']         = xsec
    sampdic['year']         = year
    sampdic['treeName']     = treeName
    sampdic['histAxisName'] = histAxisName
    sampdic['options']      = options
    sampdic['era']          = era

    # load files and get list of wc names from first file
    # files = utils.get_files(path, recursive=True)
    all_files = []
    for p in paths:
        found_files = utils.get_files(p, recursive=True)
        all_files.extend(found_files)

    if not all_files:
        print(f"No files found in the provided paths: {paths}")
        return

    # wc_names = utils.get_list_of_wc_names(files[0])
    wc_names = utils.get_list_of_wc_names(all_files[0])


    # initialize empty counters
    nevents = 0
    n_gen_events = 0
    n_sum_of_weights = 0
    is_data_lst = []

    # loop through files to get nevents and sow, check is_data
    # for f in files: 
    for f in all_files:
        # i_events, i_gen_events, i_sum_of_weights, is_data = utils.get_info(f, treeName)
        i_events, i_gen_events, i_sum_of_weights, sow_lhe_wgts, is_data = utils.get_info(f, treeName)
        nevents += i_events
        n_gen_events += i_gen_events
        n_sum_of_weights += i_sum_of_weights
        is_data_lst.append(is_data)

    if len(set(is_data_lst)) != 1:
        raise Exception("ERROR: There are a mix of files that are data and mc")
    else:
        is_data = is_data_lst[0]

    sampdic['WCnames'] = wc_names
    sampdic['files'] = [f"/store{fname.split("store")[1]}" for fname in all_files]
    sampdic['nEvents'] = nevents
    sampdic['nSumOfWeights'] = n_sum_of_weights
    sampdic['isData'] = is_data
    # sampdic['path'] = f"/store{path.split("store")[1]}"
    sampdic['path'] = [f"/store{p.split('store')[1]}" if "store" in p else p for p in paths]

    if not outname.endswith('.json'): 
        outname += '.json'

    # save sampdic in json file
    with open(outputFile, "w") as outfile:
        json.dump(sampdic, outfile, indent=4)
        print(f"\n New json file: {outputFile}") 

    print("If this is an EFT sample, the sow is incorrect and needs to be replaced with the sum of event weights at the SM")
    print("If this is a skimmed sample, REPLACE the nEvents and nSumOfWeights with the unskimmed number!")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create json file with list of samples and metadata')
    # parser.add_argument('path'              , default=''           , help = 'Path to local directory')
    parser.add_argument('paths', nargs='+', help = 'Space-separated list of paths to local directories')
    parser.add_argument('--xsec','-x'       , default=1, type=float, help = 'Cross section (number or file to read)')
    parser.add_argument('--year','-y'       , default=-1           , help = 'Year')
    parser.add_argument('--treename'        , default='Events'     , help = 'Name of the tree')
    parser.add_argument('--histAxisName'    , default=''           , help = 'Name for the samples axis of the coffea hist', required=True)
    parser.add_argument('--era'             , default=''           , help = 'Era Name') #Needed for Era dependency in Run3
    parser.add_argument('--outname','-o'    , default='temp'       , help = 'Out name of the json file', required=True)
    parser.add_argument('--options'         , default=''           , help = 'Sample-dependent options to pass to your analysis')

    args = parser.parse_args()
    # path         = args.path
    paths        = args.paths
    xsec         = args.xsec
    year         = args.year
    treeName     = args.treename
    era          = args.era
    histAxisName = args.histAxisName
    outname      = args.outname
    options      = args.options

    make_sample_json(xsec, year, treeName, histAxisName, options, era, paths, outname)
