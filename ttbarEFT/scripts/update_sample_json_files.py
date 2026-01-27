import os
import json
import argparse
import topcoffea.modules.utils as utils

def update_json(path, refpath, outname):

    # check that output json file doesn't already exist
    outputFile = outname+'.json'
    assert not os.path.exists(outputFile), "Output file will not be overwritten! Backup the file and try again."
    
    with open(refpath, 'r') as file:
        refjson = json.load(file)

    # save input arguments
    sampdic = {}
    sampdic['xsec']         = refjson['xsec']
    sampdic['year']         = refjson['year']
    sampdic['treeName']     = refjson['treeName']
    sampdic['histAxisName'] = refjson['histAxisName']
    sampdic['options']      = refjson['options']

    newjson_dict = {}
    # load files and get list of wc names from first file
    files = utils.get_files(path)
    wc_names = utils.get_list_of_wc_names(files[0])
    treeName = "Events"

    # initialize empty counters
    nevents = 0
    n_gen_events = 0
    n_sum_of_weights = 0
    is_data_lst = []

    # loop through files to get nevents and sow, check is_data
    for f in files: 
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

    newjson_dict['WCnames'] = wc_names
    newjson_dict['isData'] = is_data

    # newjson_dict['nEvents'] = nevents
    # newjson_dict['nSumOfWeights'] = n_sum_of_weights

    newjson_dict['files'] = files
    newjson_dict['path'] = path

    for key in ['WCnames', 'isData']:
        assert refjson[key] == newjson_dict[key], f"{key} entry is not the same in the old and new jsons."

    sampdic['WCnames']      = refjson['WCnames']
    # sampdic['files']        = newjson_dict['files']
    sampdic['files']        = [f"/store{fname.split("store")[1]}" for fname in newjson_dict['files']]
    sampdic['nEvents']      = refjson['nEvents']
    sampdic['nSumOfWeights']= refjson['nSumOfWeights']
    sampdic['isData']       = refjson['isData']
    # sampdic['path']         = refjson['path']

    try:
        sampdic['path']     = refjson['path']
    except: 
        sampdic['path']     = os.path.split(refjson['files'][0])[0]

    if not outname.endswith('.json'): 
        outname += '.json'

    # save sampdic in json file
    with open(outputFile, "w") as outfile:
        json.dump(sampdic, outfile, indent=4)
        print(f"\n New json file: {outputFile}") 

    print("If this is an EFT sample, the sow is incorrect and needs to be replaced with the sum of event weights at the SM")
    print("If this is a skimmed sample, check that the nEvents and nSumOfWeights matches the original unskimmed number! \n\n")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create json file with list of samples and metadata')
    parser.add_argument('--path'            , default=''           , help = 'Path to local directory')
    parser.add_argument('--refjson'         , required=True        , help = 'Path to reference json')
    parser.add_argument('--outname','-o'    , required=True        , help = 'Out name of the json file')

    args = parser.parse_args()
    path         = args.path
    refpath      = args.refjson
    outname      = args.outname

    update_json(path, refpath, outname)
