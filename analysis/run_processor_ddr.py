#!/usr/bin/env python
import argparse
import json
import time
import cloudpickle
import gzip
import os
import importlib

import numpy as np
from coffea import dataset_tools
from coffea import processor
from coffea.nanoevents import NanoAODSchema

from topcoffea.modules import utils
import topcoffea.modules.remote_environment as remote_environment

from ddr import CoffeaDynamicDataReduction
import ndcctools.taskvine as vine

LST_OF_KNOWN_EXECUTORS = ["taskvine", "iterative"]
results_dir = "/users/hnelson2/ddr_testing/"

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('jsonFiles'        , nargs='?', help = 'Json file(s) containing files and metadata')
    # parser.add_argument('--executor','-x'  , default='work_queue', help = 'Which executor to use')
    parser.add_argument('--prefix', '-r'   , nargs='?', default='', help = 'Prefix or redirector to look for the files')
    parser.add_argument('--pretend'        , action='store_true', help = 'Read json files but, not execute the analysis')
    parser.add_argument('--nworkers','-n' , default=8  , help = 'Number of workers')
    parser.add_argument('--chunksize','-s', default=100000  , help = 'Number of events per chunk')
    parser.add_argument('--nchunks','-c'  , default=None  , help = 'You can choose to run only a number of chunks')
    parser.add_argument('--outname','-o'  , default='histos', help = 'Name of the output file with histograms')
    parser.add_argument('--treename'      , default='Events', help = 'Name of the tree inside the files')
    parser.add_argument('--wc-list', action='extend', nargs='+', help = 'Specify a list of Wilson coefficients to use in filling histograms.')
    parser.add_argument('--hist-list', action='extend', nargs='+', help = 'Specify a list of histograms to fill.')
    parser.add_argument('--port', default='9123-9130', help = 'Specify the Work Queue port. An integer PORT or an integer range PORT_MIN-PORT_MAX.')
    parser.add_argument('--processor', '-p', default='analysis_processor.py', help='Specify processor file name')

    args        = parser.parse_args()
    jsonFiles   = args.jsonFiles
    # executor    = args.executor
    prefix      = args.prefix
    pretend     = args.pretend
    nworkers    = int(args.nworkers)
    chunksize   = int(args.chunksize)
    nchunks     = int(args.nchunks) if not args.nchunks is None else args.nchunks
    outname     = args.outname
    treename    = args.treename
    wc_lst      = args.wc_list if args.wc_list is not None else []
    proc_file   = args.processor
    proc_name   = args.processor[:-3]
    hist_lst    = args.hist_list

    print("\n running with processor: ", proc_file, '\n')
    print("\n jsonFiles argument: ", jsonFiles, '\n')

    analysis_processor = importlib.import_module(proc_name)

    # Check if we have valid options
    # if executor not in LST_OF_KNOWN_EXECUTORS:
    #     raise Exception(f"The \"{executor}\" executor is not known. Please specify an executor from the known executors ({LST_OF_KNOWN_EXECUTORS}). Exiting.")

    # if executor in ["taskvine"]:
    # construct wq port range
    port = list(map(int, args.port.split('-')))
    if len(port) < 1:
        raise ValueError("At least one port value should be specified.")
    if len(port) > 2:
        raise ValueError("More than one port range was specified.")
    if len(port) == 1:
        # convert single values into a range of one element
        port.append(port[0])


    # Load samples from json and setup the inputs to the processor
    samplesdict = {}
    allInputFiles = []

    def LoadJsonToSampleName(jsonFile, prefix):
        sampleName = jsonFile if not '/' in jsonFile else jsonFile[jsonFile.rfind('/')+1:]
        if sampleName.endswith('.json'): sampleName = sampleName[:-5]
        with open(jsonFile) as jf:
            samplesdict[sampleName] = json.load(jf)
            samplesdict[sampleName]['redirector'] = prefix

    # print(f"isinstance(jsonFiles, str): {isinstance(jsonFiles, str)}")

    if isinstance(jsonFiles, str) and ',' in jsonFiles:
        jsonFiles = jsonFiles.replace(' ', '').split(',')
    elif isinstance(jsonFiles, str):
        jsonFiles = [jsonFiles]

    # print(f"jsonFiles: {jsonFiles}")
    
    for jsonFile in jsonFiles:
        if os.path.isdir(jsonFile):
            if not jsonFile.endswith('/'): jsonFile+='/'
            for f in os.path.listdir(jsonFile):
                if f.endswith('.json'): allInputFiles.append(jsonFile+f)
        else:
            allInputFiles.append(jsonFile)

    # print(f"allInputFiles: {allInputFiles}")

    # Read from cfg files
    for f in allInputFiles:
        if not os.path.isfile(f):
            raise Exception(f'[ERROR] Input file {f} not found!')
        # This input file is a json file, not a cfg
        if f.endswith('.json'):
            LoadJsonToSampleName(f, prefix)
        # Open cfg files
        else:
            with open(f) as fin:
                print(' >> Reading json from cfg file...')
                lines = fin.readlines()
                for l in lines:
                    if '#' in l:
                        l=l[:l.find('#')]
                    l = l.replace(' ', '').replace('\n', '')
                    if l == '': continue
                    if ',' in l:
                        l = l.split(',')
                        for nl in l:
                            if not os.path.isfile(l):
                                prefix = nl
                            else:
                                LoadJsonToSampleName(nl, prefix)
                    else:
                        if not os.path.isfile(l):
                            prefix = l
                        else:
                            LoadJsonToSampleName(l, prefix)

    ### NEW PREPROCESSING ### 
    data_dict = {}
    for sname in samplesdict.keys():

        files_dict = {}
        for f in samplesdict[sname]['files']:
            files_dict[f] = {'object_path': 'Events'}
        
        data_dict[sname] = {'files': files_dict}

    ### END OF NEW PREPROCESSING ###

    # Extract the list of all WCs, as long as we haven't already specified one.
    if len(wc_lst) == 0:
        for k in samplesdict.keys():
            for wc in samplesdict[k]['WCnames']:
                if wc not in wc_lst:
                    wc_lst.append(wc)

    if len(wc_lst) > 0:
        # Yes, why not have the output be in correct English?
        if len(wc_lst) == 1:
            wc_print = wc_lst[0]
        elif len(wc_lst) == 2:
            wc_print = wc_lst[0] + ' and ' + wc_lst[1]
        else:
            wc_print = ', '.join(wc_lst[:-1]) + ', and ' + wc_lst[-1]
            print('Wilson Coefficients: {}.'.format(wc_print))
    
    else:
        print('No Wilson coefficients specified')

    mgr = vine.Manager(
        port=port,
        name= f"{os.environ['USER']}-ddr-coffea",
        # environment_file='/users/hnelson2/ttbarEFT-coffea2025/analysis/topeft-envs/env_spec_03e46e41_edit_HEAD.tar.gz',
        # filepath= f"/tmp/{os.environ['USER']}/vine-tmp/",
        # extra_input_files=[proc_file],
    )

    mgr.tune("hungry-minimum", 1)
    mgr.enable_monitoring(watchdog=False)

    # Check if the X509 proxy file exists
    x509_proxy = f"/tmp/x509up_u{os.getuid()}"
    if not os.path.exists(x509_proxy):
        print(
            f"Warning: X509 proxy file {args.x509_proxy} does not exist. Setting to None."
        )
        x509_proxy = None

    ddr = CoffeaDynamicDataReduction(
        mgr, #taskvine manager,
        # data={'test_dataset': {
        #     "files": {
        #         "/cms/cephfs/data/store/user/hnelson2/mc/central_ttbar_nanoAOD/UL17/nanoAOD_TTto2L2Nu_1Jets_smeft_MTT_0to700/NAOD-00000_2262.root": 
        #             {'object_path': "Events",
        #              'num_entries': 1000},
        #         }
        #     }
        # },
        data = dataset_tools.preprocess(data_dict)[0], #only grab out_available not out_updated 
        processors = {
            "ee_chan": analysis_processor.AnalysisProcessor(samplesdict,'ee', wc_lst, hist_lst),
            "mm_chan": analysis_processor.AnalysisProcessor(samplesdict,'mm', wc_lst, hist_lst),
        },
        schema=NanoAODSchema,
        remote_executor_args={"scheduler": "threads"},
        checkpoint_accumulations=False,
        resources_processing={"cores": 1},
        resources_accumualting={"cores": 2},
        results_directory=results_dir,
        x509_proxy=x509_proxy,
    )

    result = ddr.compute()
    pprint.pprint(result)
