#!/usr/bin/env python
import argparse
import json
import yaml
import time
import cloudpickle
import gzip
import os
import importlib
import pprint

import numpy as np
from coffea import processor
from coffea.nanoevents import NanoAODSchema


from coffea.processor.test_items import NanoEventsProcessor
from dynamic_data_reduction import preprocess, CoffeaDynamicDataReduction

from topcoffea.modules import utils
import topcoffea.modules.remote_environment as remote_environment

import ndcctools.taskvine as vine


def get_filename_from_path(filename):
    full_file_name = os.path.basename(filename)
    base, extension = os.path.splitext(full_file_name)
    return base


def load_json_to_sampledict(inputFile, prefix):

    temp_dict = {}
    json_dict = read_json_file(inputFile)
    sample_name = get_filename_from_path(inputFile)
    temp_dict[sample_name] = json_dict
    temp_dict[sample_name]['redirector'] = prefix

    return temp_dict


def preprocessing_for_taskvine(samplesdict):
    flist = {}
    for sname in samplesdict.keys():
        redirector = samplesdict[sname]['redirector']
        flist[sname] = [(redirector+f) for f in samplesdict[sname]['files']]

        return flist


def preprocessing_for_ddr(sample):
    
    files_dict = {}
    for f in sample['files']:
        fname = sample['redirector']+f
        files_dict[fname] = {'object_path': 'Events'}
    
    return {'files': files_dict}


def read_json_file(filename):
    with open(filename) as f:
        return json.load(f)


def read_yaml_file(filename): 
    with open(filename) as f: 
        return yaml.safe_load(f)

if __name__ == '__main__':

    # TODO: make this an input argument with a default or make it based on --outname
    results_dir = "/users/hnelson2/ddr_testing/"

    #TODO: add IterativeExecutor Options
    known_executors = ['iterative', 'ddr']

    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('inputFile'        , nargs='?', help = 'Json file(s) containing files and metadata')
    parser.add_argument('--executor','-x'  , default='work_queue', help = 'Which executor to use')
    parser.add_argument('--prefix', '-r'   , nargs='?', default='', help = 'Prefix or redirector to look for the files')
    # parser.add_argument('--pretend'        , action='store_true', help = 'Read json files but, not execute the analysis')
    # parser.add_argument('--nworkers','-n' , default=8  , help = 'Number of workers')
    parser.add_argument('--chunksize','-s', default=100000  , help = 'Number of events per chunk')
    parser.add_argument('--nchunks','-c'  , default=None  , help = 'You can choose to run only a number of chunks')
    parser.add_argument('--outname','-o'  , default='histos', help = 'Name of the output file with histograms')
    parser.add_argument('--treename'      , default='Events', help = 'Name of the tree inside the files')
    parser.add_argument('--wc-list', action='extend', nargs='+', help = 'Specify a list of Wilson coefficients to use in filling histograms.')
    parser.add_argument('--hist-list', action='extend', nargs='+', help = 'Specify a list of histograms to fill.')
    parser.add_argument('--port', default='9123-9130', help = 'Specify the Work Queue port. An integer PORT or an integer range PORT_MIN-PORT_MAX.')
    parser.add_argument('--processor', '-p', default='analysis_processor.py', help='Specify processor file name')

    args        = parser.parse_args()
    inputFile   = args.inputFile
    executor    = args.executor
    prefix      = args.prefix
    # pretend     = args.pretend
    # nworkers    = int(args.nworkers)
    chunksize   = int(args.chunksize)
    nchunks     = int(args.nchunks) if not args.nchunks is None else args.nchunks
    outname     = args.outname
    treename    = args.treename
    wc_lst      = args.wc_list if args.wc_list is not None else []
    proc_file   = args.processor
    proc_name   = args.processor[:-3]
    hist_lst    = args.hist_list
    ports       = args.port

    print("\n running with processor: ", proc_file, '\n')

    analysis_processor = importlib.import_module(proc_name)

    # Check if input is json or yaml
    if inputFile.endswith('.json'):
        isJson = True 
        isYaml = False
    elif (inputFile.endswith('.yaml')) or (inputFile.endswith('.yml')):
        isJson = False
        isYaml = True
    else:
        raise ValueError(f"Expects a .json, .yaml, or .yml for the input file. inputFile ={inputFile}")

    # Check if we have valid options
    if executor not in known_executors:
        raise Exception(f"The \"{executor}\" executor is not known. Please specify an executor from the known executors ({known_executors}). Exiting.")

    if executor in ["ddr"]:
        # construct wq port range
        port = list(map(int, ports.split('-')))
        if len(port) < 1:
            raise ValueError("At least one port value should be specified.")
        if len(port) > 2:
            raise ValueError("More than one port range was specified.")
        if len(port) == 1:
            # convert single values into a range of one element
            port.append(port[0])

        print(f"port: {port}")

    ### Fill Samples Dictionary ### 
    samplesdict = {}

    if isJson: 
        samplesdict.update(load_json_to_sampledict(inputFile, prefix))

    elif isYaml: 
        yaml_dict = read_yaml_file(inputFile)
        redirector = yaml_dict['redirector']
        jsonFiles = yaml_dict['jsonFiles']

        print(f"\n\n redirector from yaml: {redirector} \n\n")

        for f in jsonFiles: 
            samplesdict.update(load_json_to_sampledict(f, redirector))

    flist = {}
    for sname in samplesdict.keys():
        #flist[sname] = samplesdict[sname]['files']
        redirector = samplesdict[sname]['redirector']
        flist[sname] = [(redirector+f) for f in samplesdict[sname]['files']]

    ### Fill WC list ### 
    if len(wc_lst) == 0:
        for k in samplesdict.keys():
            for wc in samplesdict[k]['WCnames']:
                if wc not in wc_lst:
                    wc_lst.append(wc)
    if len(wc_lst) > 0: 
        print(f"Wilson Coefficients: {wc_lst}")
    else: 
        print(f"Wilson Coefficients: NONE SPECIFIED")

    # Run the processor and get the output
    tstart = time.time()

    if executor == 'ddr': 

        mgr = vine.Manager(
            port=port, 
            name=f"{os.environ['USER']}-ddr-coffea",
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

        # print("samplesdict:")
        # pprint.pprint(samplesdict)
        # print("\n\n")

        data = {}
        for sname in samplesdict.keys():
            data[sname] = preprocessing_for_ddr(samplesdict[sname])
        # print("data for ddr:")
        # pprint.pprint(data)
        # print("\n\n")

        print("Original data spec:")
        for dataset_name, dataset_info in data.items():
            print(f"  {dataset_name}:")
            for file_path, file_info in dataset_info["files"].items():
                print(f"    {file_path}: {file_info}")

        print("\nPreprocessing data with TaskVine...")
        preprocessed_data = preprocess(
            manager=mgr,
            data=data,
            tree_name="Events",
            timeout=30,
            max_retries=3,
            show_progress=True,
            batch_size=5,
        )

        with open("{inputFile}_preprocessed.json", "w") as f:
            json.dump(preprocessed_data, f, indent=2)

        print("\nPreprocessed data spec:")
        for dataset_name, dataset_info in preprocessed_data.items():
            print(f"  {dataset_name}:")
            for file_path, file_info in dataset_info["files"].items():
                print(f"    {file_path}: {file_info}")

        ### Dynamic Data Reduction ### 
        ddr = CoffeaDynamicDataReduction(
            mgr, #taskvine manager,
            data = preprocessed_data,
            processors = {
                "ee_chan": analysis_processor.AnalysisProcessor(samplesdict,'ee', wc_lst, hist_lst),
                "mm_chan": analysis_processor.AnalysisProcessor(samplesdict,'mm', wc_lst, hist_lst),
            },
            extra_files = [proc_file],
            schema=NanoAODSchema,
            accumulator=analysis_processor.AnalysisProcessor,
            resources_processing={"cores": 1},
            resources_accumualting={"cores": 2},
            results_directory=results_dir,
            verbose=True,
            x509_proxy=x509_proxy,
        )

        hists = ddr.compute()

        print(f"\n\n resulting hists: ")
        pprint.pprint(hists)


    elif executor == 'iterative': 

        flist = preprocessing_for_taskvine(samplesdict)
        proc_instance = analysis_processor.AnalysisProcessor(samplesdict, 'ee', wc_lst,hist_lst)
        exec_instance = processor.IterativeExecutor()
        runner = processor.Runner(exec_instance, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks)
        hists = runner(fileset=flist, processor_instance=proc_instance, treename=treename)


    ### Save Output ###
    outpath = '.'
    out_pkl_file = os.path.join(outpath, f"{outname}.pkl.gz")
    print(f"\n\n Saving output in {out_pkl_file}")
    with gzip.open(out_pkl_file, "wb") as fout:
        cloudpickle.dump(hists, fout)
        print(f"Done! \n\n")


    tend = time.time()
    print(f"\n\n Total processing time: {tend-tstart}")
