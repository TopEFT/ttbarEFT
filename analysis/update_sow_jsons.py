# /usr/bin/env python
import time
import argparse
import gzip
import pickle
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
from topcoffea.modules.update_json import update_json
import topcoffea.modules.remote_environment as remote_environment

import ndcctools.taskvine as vine


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
    known_executors = ['iterative', 'ddr', 'test']

    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('inputFile'        , nargs='?', help = 'Json file(s) containing files and metadata')
    parser.add_argument('--executor','-x'  , default='work_queue', help = 'Which executor to use')
    parser.add_argument('--prefix', '-r'   , nargs='?', default='file:///cms/cephfs/data/', help = 'Prefix or redirector to look for the files')
    # parser.add_argument('--pretend'        , action='store_true', help = 'Read json files but, not execute the analysis')
    # parser.add_argument('--nworkers','-n' , default=8  , help = 'Number of worvi kers')
    parser.add_argument('--chunksize','-s', default=100000  , help = 'Number of events per chunk')
    parser.add_argument('--nchunks','-c'  , default=None  , help = 'You can choose to run only a number of chunks')
    parser.add_argument('--outname','-o'  , default='histos', help = 'Name of the output file with histograms')
    parser.add_argument('--treename'      , default='Events', help = 'Name of the tree inside the files')
    parser.add_argument('--wc-list', action='extend', nargs='+', help = 'Specify a list of Wilson coefficients to use in filling histograms.')
    parser.add_argument('--hist-list', action='extend', nargs='+', help = 'Specify a list of histograms to fill.')
    parser.add_argument('--port', default='9123-9130', help = 'Specify the Work Queue port. An integer PORT or an integer range PORT_MIN-PORT_MAX.')
    parser.add_argument('--processor', '-p', default='analysis_processor.py', help='Specify processor file name')
    parser.add_argument('--relpath', default='../input_samples/', help='relative path to jsons (equiv to relative path in yaml)')

    args        = parser.parse_args()
    inputFile   = args.inputFile
    executor    = args.executor
    prefix      = args.prefix
    chunksize   = int(args.chunksize)
    nchunks     = int(args.nchunks) if not args.nchunks is None else args.nchunks
    outname     = args.outname
    treename    = args.treename
    wc_lst      = args.wc_list if args.wc_list is not None else []
    proc_file   = args.processor
    proc_name   = args.processor[:-3]
    hist_lst    = args.hist_list
    ports       = args.port
    relpath     = args.relpath


    ##############################
    ###### Initial Settings ######
    ##############################

    tstart = time.time()

    print("\n\nrunning with processor: ", proc_file, '\n')
    analysis_processor = importlib.import_module(proc_name)

    ### Check json or yaml ###
    if inputFile.endswith('.json'):
        isJson, isYaml = True, False
    elif (inputFile.endswith('.yaml')) or (inputFile.endswith('.yml')):
        isJson, isYaml = False, True
    else:
        raise ValueError(f"Unsupported file format: {inputFile} Use a .json, .yaml, or .yml")

    # Check if we have valid options
    if executor not in known_executors:
        raise Exception(f"The \"{executor}\" executor is not known. Please specify an executor from the known executors ({known_executors}). Exiting.")


    ###############################
    ###### Create Input Data ###### 
    ###############################

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

    ### Create Input Data dictionary for (pre)processing ###
    input_data = data_for_preprocessing(samplesdict)

    print(f"samplesdict: {samplesdict}")

    ##################################
    ###### Run Using VineReduce ###### 
    ##################################

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

        # env = remote_environment.get_environment(
        #     extra_pip_local = {"ttbarEFT": ["ttbarEFT", "setup.py"],
        #                         "dynamic_data_reduction": []},
        # )
        # sched_for_pre = partial(mgr, environment=env)

        # get X509 proxy file
        x509_proxy = f"/tmp/x509up_u{os.getuid()}"
        if not os.path.exists(x509_proxy):
            print(f"Warning: X509 proxy file {args.x509_proxy} does not exist. Setting to None.")
            x509_proxy = None
        else: 
            shutil.copy(x509_proxy, "./proxy.pem")

        # create preprocessed dataset 
        preprocessed_data_path = f"{os.path.splitext(inputFile)[0]}_preprocessed.json"  # name for preprocessed data

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
            # save_to_file = inputFile.removesuffix(".json").removesuffix(".yml").removesuffix(".yaml"), #only works python3.9 and above
            save_to_file = preprocessed_data_path,
        )

        ### Dynamic Data Reduction ### 
        print(f"\n\nProcessing data with VineReduce...")
        ddr = CoffeaDynamicDataReduction(
            mgr, #taskvine manager,
            data = preprocessed_data,
            processors = {
                "sow": analysis_processor.AnalysisProcessor(samples=samplesdict, wc_names_lst=wc_lst, hist_lst=hist_lst),
            },
            extra_files = [proc_file, "/users/hnelson2/ttbarEFT-coffea2025/ttbarEFT/params/channels.json", "proxy.pem"],
            schema=NanoAODSchema,
            max_task_retries= 20, # default=10
            step_size=800000, #equivalent to chunksize, default=100k
            resources_processing={"cores": 1},
            resources_accumulating={"cores": 1},
            results_directory=results_dir,
            verbose=True,
            x509_proxy=x509_proxy,
        )
        ddr.environment_variables["X509_USER_PROXY"] = "proxy.pem"
       
        # The output of ddr_hists is separate for each channel 
        ddr_hists = ddr.compute()

        hists = {}
        # loop through the channels (items in processors from ddr)
        for ch in ddr_hists.keys():
            hists[ch] = {}
            variables = list(ddr_hists[ch][next(iter(ddr_hists[ch]))].keys())   # get list of variables present in this channel

            for var in variables: 
                h_list = [ddr_hists[ch][dataset][var] for dataset in ddr_hists[ch]]
                hists[ch][var] = functools.reduce(operator.add, h_list)         # combines different processes into the same histEFT for the same variable

        print(f"Computing done!")


    ##########################################
    ###### Run Using Iterative Executor ###### 
    ##########################################

    elif executor == 'iterative': 

        proc_instance = analysis_processor.AnalysisProcessor(samples=samplesdict, wc_names_lst=wc_lst, hist_lst=hist_lst)
        exec_instance = processor.IterativeExecutor()
        runner = processor.Runner(exec_instance, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks)
        hists = runner(fileset=input_data, processor_instance=proc_instance, treename=treename)


    ##########################
    ###### Save Outputs ###### 
    ##########################

    outpath = '.'
    out_pkl_file = os.path.join(outpath, f"{outname}.pkl.gz")
    print(f"\n\n Saving output to {out_pkl_file}")
    with gzip.open(out_pkl_file, "wb") as fout:
        cloudpickle.dump(hists, fout)
        print(f"Done! \n\n")


    ###############################
    ###### Update JSON files ###### 
    ###############################

    if executor == 'ddr':

        WEIGHTS_NAME_LST = [
            # 'nom',
            'ISRUp', 'ISRDown',
            'FSRUp', 'FSRDown',
            'renormUp', 'renormDown',
            'factUp', 'factDown',
            'renormfactUp', 'renormfactDown',
        ]

        variables = hists['sow'].keys()
        processes = list(hists['sow']['SumOfWeights'].axes['process'])

        for proc in processes: 
            json_path = os.path.join(relpath, f"{proc}.json")

            updates = {}
            for weight in WEIGHTS_NAME_LST:
                updates[f"nSumOfWeights_{weight}"] = float(hists['sow'][f"SumOfWeights_{weight}"][{'process':proc}].as_hist({}).values()[0])

            # for just ttbar samples save the top pt sow
            if (('TT01j2l' in proc) or ('TTTo2L2Nu' in proc)) and ('SumOfWeights_toppt' in variables): 
                updates[f"nSumOfWeights_toppt"] = float(hists['sow'][f"SumOfWeights_toppt"][{'process':proc}].as_hist({}).values()[0])

            # add PDF sow variations
            if 'sow_LHEPDFweights' in variables:
                PDFweights = hists['sow']['sow_LHEPDFweights'][{'process':proc}]
                updates = {}
                for i in PDFweights.axes['PDFindex']:
                    sow = PDFweights[{'PDFindex':i}].values()[0]
                    updates[f"nSumOfWeights_LHEPDFweights{i}"] = sow

            update_json(json_path,dry_run=False,verbose=True, **updates)

            try: 
                print(f"sow at SM: {hists['sow']['SumOfWeights'][proc].as_hist({}).values()[0]}")
            except: 
                print(f"error in sow print statement")

        print(f"\n\n WARNING: This script does not update the base sum of weights. If this is an EFT sample, update the sow manually using the pkl this script produces")


    tend = time.time()
    print(f"\n\nTotal processing time: {tend-tstart}")
