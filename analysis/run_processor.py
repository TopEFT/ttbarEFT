#!/usr/bin/env python
import argparse
import json
import time
import cloudpickle
import gzip
import os
import cProfile
import importlib

import numpy as np
from coffea import processor
from coffea.nanoevents import NanoAODSchema

from topcoffea.modules import utils
import topcoffea.modules.remote_environment as remote_environment

LST_OF_KNOWN_EXECUTORS = ["futures", "work_queue"]
proc_options = ["net_disc", "kinematics_processor", "weight_processor", "analysis_processor"]

def main():
    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('jsonFiles'        , nargs='?', help = 'Json file(s) containing files and metadata')
    parser.add_argument('--executor','-x'  , default='work_queue', help = 'Which executor to use')
    parser.add_argument('--prefix', '-r'   , nargs='?', default='', help = 'Prefix or redirector to look for the files')
    parser.add_argument('--nworkers','-n' , default=8  , help = 'Number of workers')
    parser.add_argument('--chunksize','-s', default=100000  , help = 'Number of events per chunk')
    parser.add_argument('--nchunks','-c'  , default=None  , help = 'You can choose to run only a number of chunks')
    parser.add_argument('--outname','-o'  , default='histos', help = 'Name of the output file with histograms')
    parser.add_argument('--treename'      , default='Events', help = 'Name of the tree inside the files')
    parser.add_argument('--hist-list', action='extend', nargs='+', help = 'Specify a list of histograms to fill.')
    parser.add_argument('--port', default='9123-9130', help = 'Specify the Work Queue port. An integer PORT or an integer range PORT_MIN-PORT_MAX.')
    parser.add_argument('--processor', '-p', default='kinematics_processor', help='Specify processor name (without .py)')
    parser.add_argument('--lnet',  '-l', default='ctq8_5.0', help='Specify linear networks as array of strings (i.e. ["ctq8_5.0", "ctq1_1.0"])')

    args        = parser.parse_args()
    jsonFiles   = args.jsonFiles
    executor    = args.executor
    prefix      = args.prefix
    nworkers    = int(args.nworkers)
    chunksize   = int(args.chunksize)
    nchunks     = int(args.nchunks) if not args.nchunks is None else args.nchunks
    outname     = args.outname
    treename    = args.treename
    proc        = args.processor
    lnet        = args.lnet

    if proc not in proc_options:
        raise Exception(f"The \"{proc}\" processor is not known. Modify the run script or choose from the list of known processors: {proc_options}. Exiting. ")

    proc_file = proc+'.py'
    print("\n running with processor: ", proc_file, '\n')

    if executor not in LST_OF_KNOWN_EXECUTORS:
        raise Exception(f"The \"{executor}\" executor is not known. Please specify an executor from the known executors ({LST_OF_KNOWN_EXECUTORS}). Exiting.")

    if executor == "work_queue":
        port = list(map(int, args.port.split('-')))
        if len(port) < 1:
            raise ValueError("At least one port value should be specified.")
        if len(port) > 2:
            raise ValueError("More than one port range was specified.")
        if len(port) == 1:
            port.append(port[0])

    if args.hist_list == ["kinematicsDNN"]:
        hist_lst = ['Lep1_pt', 'Lep2_pt', 'Lep1_eta',
                    'Lep2_eta', 'Lep1_phi', 'Lep2_phi', 
                    'nJet30', 'jet1_pt', 'jet2_pt']
    elif args.hist_list == ["net_disc"]:
        hist_lst = ['net_disc']
    elif args.hist_list == ["kinematics"]:
        hist_lst = ["tops_pt", "avg_top_pt", "l0pt", 
                    "dr_leps", "ht", "jets_pt", 
                    "j0pt", "njets", "mtt", "mll"]
    else:
        hist_lst = args.hist_list

    samplesdict = {}
    allInputFiles = []

    def LoadJsonToSampleName(jsonFile, prefix):
        sampleName = jsonFile if not '/' in jsonFile else jsonFile[jsonFile.rfind('/')+1:]
        if sampleName.endswith('.json'): sampleName = sampleName[:-5]
        with open(jsonFile) as jf:
            samplesdict[sampleName] = json.load(jf)
            samplesdict[sampleName]['redirector'] = prefix

    if isinstance(jsonFiles, str) and ',' in jsonFiles:
        jsonFiles = jsonFiles.replace(' ', '').split(',')
    elif isinstance(jsonFiles, str):
        jsonFiles = [jsonFiles]
    for jsonFile in jsonFiles:
        if os.path.isdir(jsonFile):
            if not jsonFile.endswith('/'): jsonFile+='/'
            for f in os.path.listdir(jsonFile):
                if f.endswith('.json'): allInputFiles.append(jsonFile+f)
        else:
            allInputFiles.append(jsonFile)

    for f in allInputFiles:
        if not os.path.isfile(f):
            raise Exception(f'[ERROR] Input file {f} not found!')
        if f.endswith('.json'):
            LoadJsonToSampleName(f, prefix)
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
        redirector = samplesdict[sname]['redirector']
        flist[sname] = [(redirector+f) for f in samplesdict[sname]['files']]
        
    wc_set = ()
    for item in samplesdict:
        for i, file_name in enumerate(samplesdict[item]['files']):
            wc_lst = utils.get_list_of_wc_names(file_name)
            if i==0:
                wc_set = set(wc_lst)
            else:
                if set(wc_lst) != wc_set:
                   raise Exception("ERROR: Not all files have same WC list")

    if executor == "work_queue":
        executor_args = {
            'master_name': '{}-workqueue-coffea'.format(os.environ['USER']),
            'port': port,
            'debug_log': 'debug.log',
            'transactions_log': 'tr.log',
            'stats_log': 'stats.log',
            'tasks_accum_log': 'tasks.log',
            'environment_file': remote_environment.get_environment(
                extra_conda=["pytorch=2.3.0", "numpy=1.23.5", "pyyaml=6.0.1"],
                extra_pip_local = {"EFTmva": ["models", "net.py"], 
                                  },
            ),
            'extra_input_files' : [proc_file],
            'retries': 5,
            'compression': 9,
            'resource_monitor': True,
            'resources_mode': 'auto',
            'chunks_per_accum': 25,
            'chunks_accum_in_mem': 2,
            'fast_terminate_workers': 0,
            'verbose': True,
            'print_stdout': False,
        }

    tstart = time.time()

    if executor == "futures":
        exec_instance = processor.FuturesExecutor(workers=nworkers)
        runner = processor.Runner(exec_instance, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks)
    elif executor ==  "work_queue":
        executor = processor.WorkQueueExecutor(**executor_args)
        runner = processor.Runner(executor, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks, skipbadfiles=False, xrootdtimeout=180)

    if proc == 'kinematics_processor':
        import kinematics_processor
        processor_instance = kinematics_processor.AnalysisProcessor(samplesdict, wc_lst, hist_lst)
    elif proc == 'weight_processor':
        import weight_processor
        processor_instance = weight_processor.AnalysisProcessor(samplesdict, wc_lst, hist_lst, lnet=lnet) 
    elif proc == 'net_disc':
        import net_disc
        processor_instance = net_disc.AnalysisProcessor(samplesdict, wc_lst, hist_lst)
    elif proc == 'analysis_processor':
        analysis_processor = importlib.import_module(proc_name)
        processor_instance = analysis_processor.AnalysisProcessor(samplesdict,wc_lst,hist_lst)

    output = runner(flist, treename, processor_instance)
    
    dt = time.time() - tstart

    if executor == "work_queue":
        print(f'Processed {nevts_total} events in {dt:.2f} seconds ({nevts_total/dt:.2f} evts/sec).')
    elif executor == "futures":
        print(f'Processing time: {dt:.2f} with {nworkers} ({dt*nworkers:.2f} cpu overall)')

    if not os.path.isdir(outname): os.system("mkdir -p %s"%outname)
    out_pkl_file = os.path.join(outname,"histeft.pkl.gz")
    print(f"\nSaving output in {out_pkl_file}...")
    with gzip.open(out_pkl_file, "wb") as fout:
        cloudpickle.dump(output, fout)
        
    
    print("Done!")

if __name__ == '__main__':
    profile = False
    if profile:
        cProfile.run('main()', filename='profile.txt')
    else:
        main()
