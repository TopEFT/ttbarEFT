# /usr/bin/env python
import time
import argparse
import gzip
import os
import shutil
import importlib
import pprint
import json
import yaml
import cloudpickle
import functools
import operator 

from coffea import processor
from coffea.nanoevents import NanoAODSchema

from coffea.processor.test_items import NanoEventsProcessor
from dynamic_data_reduction import preprocess, CoffeaDynamicDataReduction

from topcoffea.modules import utils
import topcoffea.modules.remote_environment as remote_environment

import ndcctools.taskvine as vine

# import reweight_processor as analysis_processor

def check_preprocessed_data(input_data, preprocessed_data_path): 
    check = False 

    if os.path.exists(preprocessed_data_path): 
        preprocessed_data = read_json_file(preprocessed_data_path)

        for sname in input_data.keys():
            required_file_list = input_data[sname]['files'].keys()

            # check that the dataset exists in the preprocessed data json
            if sname in preprocessed_data.keys():
                preprocessed_file_list = preprocessed_data[sname]['files'].keys()
                check = (sorted(required_file_list) == sorted(preprocessed_file_list))

            # if dataset doesn't exist, the jsons are not the same and the check fails
            else: 
                check = False

    # if the preprocessed data json doesn't exist, check fails
    else: 
        check = False

    return check


def get_filename_from_path(filename):
    full_file_name = os.path.basename(filename)
    base, extension = os.path.splitext(full_file_name)

    return base


def load_json_to_samplesdict(inputFile, prefix):
    samplesdict = {}
    json_dict = read_json_file(inputFile)
    sample_name = get_filename_from_path(inputFile)
    samplesdict[sample_name] = json_dict
    samplesdict[sample_name]['redirector'] = prefix

    return samplesdict


def preprocessing_for_taskvine(samplesdict):
    flist = {}
    for sname in samplesdict.keys():
        redirector = samplesdict[sname]['redirector']
        flist[sname] = [(redirector+f) for f in samplesdict[sname]['files']]

    return flist


def data_for_preprocessing(samplesdict):
    data = {}
    for sname, sample in samplesdict.items():
        files_dict = {}
        metadata = dict(sample)
        del metadata['files']
        for f in sample['files']:
            fname = sample['redirector']+f
            files_dict[fname] = {'object_path': 'Events'}

        data[sname] = {'files': files_dict, 'metadata': metadata}
    
    return data


def read_json_file(filename):
    with open(filename) as f:
        return json.load(f)


def read_yaml_file(filename): 
    with open(filename) as f: 
        return yaml.safe_load(f)


if __name__ == '__main__':

    # TODO: make this an input argument with a default or make it based on --outname
    results_dir = f"/users/{os.environ['USER']}/ddr_coffea_test/"

    timestamp = time.strftime('%Y%m%d_%H%M', time.localtime())

    #TODO: add IterativeExecutor Options
    known_executors = ['iterative', 'ddr']

    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('inputFile',            nargs='?', help = 'Json file(s) containing files and metadata')
    parser.add_argument('--processor', '-p',    default='analysis_processor.py', help='Specify processor file name')
    parser.add_argument('--executor','-x',      default='work_queue', help = 'Which executor to use')
    parser.add_argument('--outname','-o',       default='histos', help = 'Name of the output file with histograms')
    parser.add_argument('--no-group',           action='store_false', default=True, dest='aggregate', help='Disable output histogram combination')
    parser.add_argument('--prefix', '-r',       nargs='?', default='', help = 'Prefix or redirector to look for the files')
    parser.add_argument('--treename',           default='Events', help = 'Name of the tree inside the files')
    parser.add_argument('--wc-list',            action='extend', nargs='+', help = 'Specify a list of Wilson coefficients to use in filling histograms.')
    parser.add_argument('--hist-list',          action='extend', nargs='+', help = 'Specify a list of histograms to fill.')
    parser.add_argument('--chunksize','-s',     default=100000  , help = 'Number of events per chunk')
    parser.add_argument('--nchunks','-c',       default=2  , help = 'Max number of chunks for IterativeExecutor')
    parser.add_argument('--port',               default='9123-9130', help = 'Specify the Work Queue port. An integer PORT or an integer range PORT_MIN-PORT_MAX.')
    
    args        = parser.parse_args()
    inputFile   = args.inputFile
    executor    = args.executor
    outname     = args.outname
    aggregate   = args.aggregate
    proc_file   = args.processor
    prefix      = args.prefix
    treename    = args.treename
    wc_lst      = args.wc_list if args.wc_list is not None else []
    hist_lst    = args.hist_list
    chunksize   = int(args.chunksize)
    nchunks     = int(args.nchunks) # if not args.nchunks is None else args.nchunks
    ports       = args.port


    print(f"running with aggregate setting: {aggregate}")
    if aggregate:
        print("True for aggregate")

    analysis_processor = importlib.import_module(proc_file[:-3])

    ### Check json or yaml ###
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

    ### Fill Samples Dictionary ### 
    samplesdict = {}

    if isJson: 
        samplesdict.update(load_json_to_samplesdict(inputFile, prefix))

    elif isYaml: 
        yaml_dict = read_yaml_file(inputFile)
        if 'jsonFiles' in yaml_dict.keys(): 
            redirector = yaml_dict['redirector']
            jsonFiles = yaml_dict['jsonFiles']

            for f in jsonFiles: 
                samplesdict.update(load_json_to_samplesdict(f, redirector))

        else: 
            for item in yaml_dict: 
                redirector = yaml_dict[item]['redirector']
                jsonFiles = yaml_dict[item]['jsonFiles']

                for f in jsonFiles: 
                    samplesdict.update(load_json_to_samplesdict(f, redirector))


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

    processor_versions = {
        "ttbar" : analysis_processor.AnalysisProcessor(samples=samplesdict),
    }

    ### RUN PROCESSOR USING VINE REDUCE ###
    if executor == 'ddr': 
        # construct port range
        port = list(map(int, ports.split('-')))
        if len(port) < 1:
            raise ValueError("At least one port value should be specified.")
        if len(port) > 2:
            raise ValueError("More than one port range was specified.")
        if len(port) == 1:
            port.append(port[0])  # convert single values into a range of one element

        # create TaskVine Manager
        mgr = vine.Manager(
            port=port, 
            name=f"{os.environ['USER']}-vineReduce-{timestamp}",
        )
        mgr.tune("hungry-minimum", 1)
        mgr.enable_monitoring(watchdog=False)
        mgr.enable_disconnect_slow_workers(5)

        # get X509 proxy file
        x509_proxy = f"/tmp/x509up_u{os.getuid()}"
        if not os.path.exists(x509_proxy):
            print(f"Warning: X509 proxy file {args.x509_proxy} does not exist. Setting to None.")
            x509_proxy = None
        else: 
            shutil.copy(x509_proxy, "./proxy.pem")

        input_data = data_for_preprocessing(samplesdict)
        preprocessed_data_path = f"{os.path.splitext(inputFile)[0]}_preprocessed.json"  # name for preprocessed data 
        # bool, checks if existing preprocessed json exists and has identical file list as input_data
        preprocessed_json_exists = check_preprocessed_data(input_data, preprocessed_data_path) 

        # use existing preprocessed json if available 
        if preprocessed_json_exists: 
            preprocessed_data = read_json_file(preprocessed_data_path)
            print(f"\n\nPreprocessed data json already exists. Using file: {preprocessed_data_path}")

        # create preprocessed dataset if not already available
        else: 
            print("\n\nPreprocessing data with TaskVine...")
            preprocessed_data = preprocess(
                manager=mgr,
                data=input_data,
                tree_name="Events",
                timeout=30,
                max_retries=5,
                show_progress=True,
                batch_size=5,
                x509_proxy=x509_proxy,
                save_to_file = preprocessed_data_path,
            )

        ### Dynamic Data Reduction ### 
        print(f"\n\nProcessing data with VineReduce...")
        ddr = CoffeaDynamicDataReduction(
            mgr, #taskvine manager,
            data = preprocessed_data,
            processors = processor_versions,
            extra_files = [proc_file, "proxy.pem"], #"/users/hnelson2/ttbarEFT-coffea2025/ttbarEFT/params/channels.json", 
            schema=NanoAODSchema,
            max_task_retries= 10, # default=10
            step_size = 10000, # 500000,
            # step_size=1000000, #equivalent to chunksize, default=100k
            resources_processing={"cores": 1, "memory":1000},
            resources_accumulating={"cores": 1},
            results_directory=results_dir,
            verbose=True,
            x509_proxy=x509_proxy,
        )
        ddr.environment_variables["X509_USER_PROXY"] = "proxy.pem"
        
        ddr_hists = ddr.compute()

        if aggregate: 
            print(f"Combining histograms from different datasets...")
            hists = {}
            
            for ch in ddr_hists.keys():
                hists[ch] = {}
                variables = list(ddr_hists[ch][next(iter(ddr_hists[ch]))].keys())

                for var in variables: 
                    h_list = [ddr_hists[ch][dataset][var] for dataset in ddr_hists[ch]]
                    hists[ch][var] = functools.reduce(operator.add, h_list)

        else: 
            hists = ddr_hists

        print(f"Computing done!")


    ### RUN PROCESSOR USING ITERATIVE EXECUTOR ###
    elif executor == 'iterative': 

        flist = preprocessing_for_taskvine(samplesdict)
        proc_instance = analysis_processor.AnalysisProcessor(samples=samplesdict)
        # proc_instance = analysis_processor.AnalysisProcessor(samples=samplesdict, lep_cat='ee', wc_names_lst=wc_lst, hist_lst=hist_lst, do_errors=do_err, syst_list=['elecID', 'muonID', 'muonISO','trigSF'])
        exec_instance = processor.IterativeExecutor()
        runner = processor.Runner(exec_instance, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks)
        hists = runner(fileset=flist, processor_instance=proc_instance, treename=treename)


    ### SAVE OUTPUT ###
    outpath = '.'
    out_pkl_file = os.path.join(outpath, f"{outname}.pkl.gz")
    print(f"\n\n Saving output to {out_pkl_file}")
    with gzip.open(out_pkl_file, "wb") as fout:
        cloudpickle.dump(hists, fout)
        print(f"Done! \n\n")

    tend = time.time()
    print(f"\n\nTotal processing time: {tend-tstart}")
