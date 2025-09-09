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

LST_OF_KNOWN_EXECUTORS = ["taskvine", "iterative"]

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('jsonFiles'        , nargs='?', help = 'Json file(s) containing files and metadata')
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
    jsonFiles   = args.jsonFiles
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
    print("\n jsonFiles argument: ", jsonFiles, '\n')

    analysis_processor = importlib.import_module(proc_name)

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


    # Load samples from json and setup the inputs to the processor
    samplesdict = {}
    allInputFiles = []

    def LoadJsonToSampleName(jsonFile, prefix):
        sampleName = jsonFile if not '/' in jsonFile else jsonFile[jsonFile.rfind('/')+1:]
        if sampleName.endswith('.json'): sampleName = sampleName[:-5]
        with open(jsonFile) as jf:
            samplesdict[sampleName] = json.load(jf)
            samplesdict[sampleName]['redirector'] = prefix

    print(f"isinstance(jsonFiles, str): {isinstance(jsonFiles, str)}")

    if isinstance(jsonFiles, str) and ',' in jsonFiles:
        jsonFiles = jsonFiles.replace(' ', '').split(',')
    elif isinstance(jsonFiles, str):
        jsonFiles = [jsonFiles]

    print(f"jsonFiles: {jsonFiles}")
    
    for jsonFile in jsonFiles:
        if os.path.isdir(jsonFile):
            if not jsonFile.endswith('/'): jsonFile+='/'
            for f in os.path.listdir(jsonFile):
                if f.endswith('.json'): allInputFiles.append(jsonFile+f)
        else:
            allInputFiles.append(jsonFile)

    print(f"allInputFiles: {allInputFiles}")

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

    flist = {}
    for sname in samplesdict.keys():
        #flist[sname] = samplesdict[sname]['files']
        redirector = samplesdict[sname]['redirector']
        flist[sname] = [(redirector+f) for f in samplesdict[sname]['files']]

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

    proc_instance = analysis_processor.AnalysisProcessor(samplesdict,wc_lst,hist_lst)

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

            'environment_file': '/users/hnelson2/ttbarEFT-coffea2025/analysis/topeft-envs/env_spec_b3473f78_edit_HEAD.tar.gz',
            # 'environment_file': remote_environment.get_environment(
            #     extra_pip_local = {"ttbarEFT": ["ttbarEFT", "setup.py"]},
            #     extra_pip=['mt2'],
            #     # extra_conda=["pytorch=2.3.1", "numpy=1.23.5"]
            # ),
            # 'extra_input_files': ["nanogen_processor.py"],
            # 'extra_input_files' : [proc_file],

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
            'print_stdout': False,
        }


    # Run the processor and get the output
    tstart = time.time()

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
    dt = time.time() - tstart

    if executor in ["taskvine"]:
        print(f"Processed {nevts_total} events in {dt} seconds ({nevts_total/dt:.2f} events/sec)")
    elif executor == "iterative":
        print(f"Processing time: {dt:.2f} seconds with {nworkers} workers ({dt*nworkers} cpu overall)")
    
    # Save the output
    outpath = "."
    if not os.path.isdir(outpath): os.system("mkdir -p %s"%outpath)
    out_pkl_file = os.path.join(outpath,outname+".pkl.gz")
    print(f"\nSaving output in {out_pkl_file}...")
    with gzip.open(out_pkl_file, "wb") as fout:
        cloudpickle.dump(output, fout)
    print("Done!")
