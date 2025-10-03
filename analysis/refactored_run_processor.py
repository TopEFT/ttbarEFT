#!/usr/bin/env python
import argparse
import json
import time
import cloudpickle
import gzip
import os
import importlib

import numpy as np
from coffea import processor
from coffea.nanoevents import NanoAODSchema

from topcoffea.modules import utils
import topcoffea.modules.remote_environment as remote_environment

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
        files_dict[f] = {'object_path': 'Events'}
    
    return {'files': files_dict}


def read_json_file(filename):
    with open(filename) as f:
        return json.load(f)


def read_yaml_file(filename): 
    with open(filename) as f: 
        return yaml.safe_load(f)

if __name__ == '__main__':

    LST_OF_KNOWN_EXECUTORS = ["taskvine", "iterative"]

    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('inputFile'        , nargs='?', help = 'Json file(s) containing files and metadata')
    parser.add_argument('--executor','-x'  , default='work_queue', help = 'Which executor to use')
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
    inputFile   = args.inputFile
    executor    = args.executor
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
    if executor not in LST_OF_KNOWN_EXECUTORS:
        raise Exception(f"The \"{executor}\" executor is not known. Please specify an executor from the known executors ({LST_OF_KNOWN_EXECUTORS}). Exiting.")

    if executor in ["taskvine"]:
        # construct wq port range
        port = list(map(int, args.port.split('-')))
        if len(port) < 1:
            raise ValueError("At least one port value should be specified.")
        if len(port) > 2:
            raise ValueError("More than one port range was specified.")
        if len(port) == 1:
            # convert single values into a range of one element
            port.append(port[0])

    ### Fill Samples Dictionary ### 
    samplesdict = {}

    if isJson: 
        samplesdict.update(load_json_to_sampledict(inputFile, prefix))

    elif isYaml: 
        yaml_dict = read_yaml_file(inputFile)
        redirector = yaml_dict['redirector']
        jsonFiles = yaml_dict['jsonFiles']

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

    if executor in ["taskvine"]:
        executor_args = {
            'manager_name': f"{os.environ['USER']}-{executor}-coffea",
            'filepath': f"/tmp/{os.environ['USER']}/vine-tmp/",

            # find a port to run work queue in this range:
            'port': port,
            # 'debug_log': 'debug.log',
            # 'transactions_log': 'tr.log',
            # 'stats_log': 'stats.log',
            # 'tasks_accum_log': 'tasks.log',

            # 'environment_file': '/users/hnelson2/ttbarEFT-coffea2025/analysis/topeft-envs/env_spec_03e46e41_edit_HEAD.tar.gz',
            # 'environment_file': remote_environment.get_environment(
            #     extra_pip_local = {"ttbarEFT": ["ttbarEFT", "setup.py"]},
            #     # extra_pip=['mt2'],
            #     # extra_conda=["pytorch=2.3.1", "numpy=1.23.5"]
            # ),
            'extra_input_files' : [proc_file],

            'retries': 10,

            # use mid-range compression for chunks results. 9 is the default for work
            # queue in coffea. Valid values are 0 (minimum compression, less memory
            # usage) to 16 (maximum compression, more memory usage).
            'compression': 8,

            # automatically find an adequate resource allocation for tasks.
            # tasks are first tried using the maximum resources seen of previously ran
            # tasks. on resource exhaustion, they are retried with the maximum resource
            # values, if specified below. if a maximum is not specified, the task waits
            # forever until a larger worker connects.
            'resource_monitor': 'measure',
            'resources_mode': 'max',

            # this resource values may be omitted when using
            # resources_mode: 'auto', but they do make the initial portion
            # of a workflow run a little bit faster.
            # Rather than using whole workers in the exploratory mode of
            # resources_mode: auto, tasks are forever limited to a maximum
            # of 8GB of mem and disk.
            #
            # NOTE: The very first tasks in the exploratory
            # mode will use the values specified here, so workers need to be at least
            # this large. If left unspecified, tasks will use whole workers in the
            # exploratory mode.
            # 'cores': 1,
            # 'disk': 8000,   #MB
            # 'memory': 10000, #MB

            # control the size of accumulation tasks. Results are
            # accumulated in groups of size chunks_per_accum, keeping at
            # most chunks_per_accum at the same time in memory per task.
            # 'chunks_per_accum': 25,
            # 'chunks_accum_in_mem': 2,

            # terminate workers on which tasks have been running longer than average.
            # This is useful for temporary conditions on worker nodes where a task will
            # be finish faster is ran in another worker.
            # the time limit is computed by multipliying the average runtime of tasks
            # by the value of 'fast_terminate_workers'.  Since some tasks can be
            # legitimately slow, no task can trigger the termination of workers twice.
            #
            # warning: small values (e.g. close to 1) may cause the workflow to misbehave,
            # as most tasks will be terminated.
            #
            # Less than 1 disables it.
            'fast_terminate_workers': 0,

            # print messages when tasks are submitted, finished, etc.,
            # together with their resource allocation and usage. If a task
            # fails, its standard output is also printed, so we can turn
            # off print_stdout for all tasks.
            'verbose': True,
            'print_stdout': True,
        }


    # Run the processor and get the output
    tstart = time.time()

    proc_dict = {
        'ee_chan': analysis_processor.AnalysisProcessor(samplesdict,'ee', wc_lst, hist_lst),
        'mm_chan': analysis_processor.AnalysisProcessor(samplesdict,'mm', wc_lst, hist_lst),
    }

    # proc_dict = {
    #     'mm_chan': analysis_processor.AnalysisProcessor(samplesdict,'mm', wc_lst, hist_lst),
    # }

    for ch, proc_instance in proc_dict.items():
        if executor == "iterative":
            exec_instance = processor.IterativeExecutor()
            runner = processor.Runner(exec_instance, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks)

        elif executor == "taskvine": 
            exec_instance = processor.TaskVineExecutor(**executor_args)
            runner = processor.Runner(
                exec_instance,
                schema=NanoAODSchema,
                chunksize=chunksize,
                maxchunks=nchunks,
                xrootdtimeout=300,
            )

        output = runner(fileset=flist, processor_instance=proc_instance, treename=treename)

        outpath = '.'
        out_pkl_file = os.path.join(outpath, f"{outname}_{ch}.pkl.gz")
        print(f"\n\n Saving output in {out_pkl_file}")
        with gzip.open(out_pkl_file, "wb") as fout:
            cloudpickle.dump(output, fout)
            print(f"Done with channel '{ch}'! \n\n")


    tend = time.time()
    print(f"\n\n Total processing time: {tend-tstart}")

