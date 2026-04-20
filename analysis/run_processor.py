#!/usr/bin/env python
import argparse, os

from coffea.nanoevents import NanoAODSchema
from dynamic_data_reduction import preprocess, CoffeaDynamicDataReduction
from importlib import import_module
from json import load
from numpy import float64, float32
from shutil import copy
from time import time
from torch.utils.data import TensorDataset
from topcoffea.modules import utils
from yaml import safe_load

import ndcctools.taskvine as vine
import topcoffea.modules.remote_environment as remote_environment

def check_preprocessed_data(input_data, preprocessed_data_path): 
    
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

def read_json_file(filename):
    with open(filename) as f:
        return load(f)

def read_yaml_file(filename): 
    with open(filename) as f: 
        return safe_load(f)
    
if __name__ == '__main__':
    known_executors = ['iterative', 'ddr']
    
    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('inputFile'        , nargs='?', help = 'Json or yaml file(s) containing files and metadata')
    parser.add_argument('--executor','-x'  , default='ddr', help = 'Which executor to use')
    parser.add_argument('--prefix', '-r'   , nargs='?', default='', help = 'Prefix or redirector to look for the files')
    parser.add_argument('--chunksize','-s', default=100000  , help = 'Number of events per chunk')
    parser.add_argument('--nchunks','-c'  , default=None  , help = 'You can choose to run only a number of chunks')
    parser.add_argument('--outname','-o'  , default='histos', help = 'Name of the output file with histograms')
    parser.add_argument('--treename'      , default='Events', help = 'Name of the tree inside the files')
    parser.add_argument('--wc-list', action='extend', nargs='+', help = 'Specify a list of Wilson coefficients to use in filling histograms.')
    parser.add_argument('--port', default='9123-9130', help = 'Specify the Work Queue port. An integer PORT or an integer range PORT_MIN-PORT_MAX.')
    parser.add_argument('--processor', '-p', default='centralGen', help='Specify processor name (without .py)')
    parser.add_argument('--dtype', '-d', default=float32, help='dtype used to build tensors')

    args        = parser.parse_args()
    inputFile   = args.inputFile
    executor    = args.executor
    prefix      = args.prefix
    chunksize   = int(args.chunksize)
    nchunks     = int(args.nchunks) if not args.nchunks is None else args.nchunks
    outname     = args.outname
    treename    = args.treename
    wc_lst      = args.wc_list if args.wc_list is not None else []
    proc        = args.processor
    ports       = args.port
    dtype       = args.dtype

    proc_file = proc+'.py'
    print("\n running with processor: ", proc_file, '\n')
    
    if proc == 'centralGen':
        import centralGen
        analysis_processor = centralGen
    elif proc == 'partons':
        import partons
        analysis_processor = partons

    isJson = inputFile.endswith('.json')
    isYaml = (inputFile.endswith('.yaml')) or (inputFile.endswith('.yml'))

    ### Check json or yaml ###
    if not isJson and not isYaml:
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
    tstart = time()

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
            name=f"{os.environ['USER']}-ddr-coffea",
        )
        mgr.tune("hungry-minimum", 1)
        mgr.enable_monitoring(watchdog=False)

        print(f'manager: {f"{os.environ['USER']}-ddr-coffea"}')
        # get X509 proxy file
        x509_proxy = f"/tmp/x509up_u{os.getuid()}"
        if not os.path.exists(x509_proxy):
            print(f"Warning: X509 proxy file {x509_proxy} does not exist. Setting to None.")
            x509_proxy = None
        else: 
            copy(x509_proxy, "./proxy.pem")

        
        # check for and create preprocessed data
        input_data = data_for_preprocessing(samplesdict)
        preprocessed_data_path = f"{os.path.splitext(inputFile)[0]}_preprocessed.json"  # name for preprocessed data 
        # bool, checks if existing preprocessed json exists and has identical file list as input_data
        preprocessed_json_exists = check_preprocessed_data(input_data, preprocessed_data_path) 

#        print(f'using {input_data}')
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
            )

        ### Dynamic Data Reduction ### 
        print(f"\n\nProcessing data with VineReduce...")
        ddr = CoffeaDynamicDataReduction(
            mgr, #taskvine manager,
            data = preprocessed_data,
            processors = {
                "tensors": analysis_processor.AnalysisProcessor(samples=samplesdict, dtype=dtype, wcs=wc_lst, outname=outname)
            },
            extra_files = [proc_file, "proxy.pem"], 
            schema=NanoAODSchema,
            max_task_retries= 20, # default=10
            step_size=800000, #equivalent to chunksize, default=100k
            resources_processing={"cores": 1},
            resources_accumulating={"cores": 1},
            results_directory=outname,
            verbose=True,
            x509_proxy=x509_proxy,
        )
        ddr.environment_variables["X509_USER_PROXY"] = "proxy.pem"
        output = ddr.compute()

        print(f"Computing done!")
        
    tend = time()
    print(f"\n\nTotal processing time: {tend-tstart}")