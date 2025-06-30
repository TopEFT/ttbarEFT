#!/usr/bin/env python
import argparse, json, os, time

from coffea import processor
from coffea.nanoevents import NanoAODSchema
from torch import float64, save

from topcoffea.modules import utils
import topcoffea.modules.remote_environment as remote_environment

def main():
    parser = argparse.ArgumentParser(description='You can customize your run')
    parser.add_argument('jsonFiles'        , nargs='?', help = 'Json file(s) containing files and metadata')
    parser.add_argument('--executor','-x'  , default='work_queue', help = 'Which executor to use')
    parser.add_argument('--prefix', '-r'   , nargs='?', default='', help = 'Prefix or redirector to look for the files')
    parser.add_argument('--nworkers','-n' , default=8  , help = 'Number of workers')
    parser.add_argument('--chunksize','-s', default=100000  , help = 'Number of events per chunk')
    parser.add_argument('--nchunks','-c'  , default=None  , help = 'You can choose to run only a number of chunks')
    parser.add_argument('--outname','-o'  , default='/scratch365/cmcgrad2/data', help = 'Name of the output file')
    parser.add_argument('--treename'      , default='Events', help = 'Name of the tree inside the files')
    parser.add_argument('--port', default='9123-9130', help = 'Specify the Work Queue port. An integer PORT or an integer range PORT_MIN-PORT_MAX.')
    parser.add_argument('--processor', '-p', default='centralGen', help='Specify processor name (without .py)')
    parser.add_argument('--dtype', '-d', default=float64, help='dtype used to build tensors')
    
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
    dtype       = args.dtype

    proc_file = proc+'.py'
    print("\n running with processor: ", proc_file, '\n')

    if executor == "work_queue":
        port = list(map(int, args.port.split('-')))
        if len(port) < 1:
            raise ValueError("At least one port value should be specified.")
        if len(port) > 2:
            raise ValueError("More than one port range was specified.")
        if len(port) == 1:
            port.append(port[0])

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
        redirector = samplesdict[sname]['redirector']
        flist[sname] = [f'{prefix}{f}' for f in samplesdict[sname]['files']]

    wc_set = ()
    for item in samplesdict:
        for i, file_name in enumerate(samplesdict[item]['files']):
            if prefix:
                wc_lst = utils.get_list_of_wc_names(f'{prefix}/{file_name}')
            else:
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
                extra_pip_local = {"ttbarEFT": ["ttbarEFT", "setup.py"]},
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
            'print_stdout': True,
        }

    tstart = time.time()

    if executor == "futures":
        exec_instance = processor.FuturesExecutor(workers=nworkers)
        runner = processor.Runner(exec_instance, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks)
    elif executor ==  "work_queue":
        executor = processor.WorkQueueExecutor(**executor_args)
        runner = processor.Runner(executor, schema=NanoAODSchema, chunksize=chunksize, maxchunks=nchunks, skipbadfiles=False, xrootdtimeout=180)

    if proc == 'ml_data':
        import ml_data
        processor_instance = ml_data.AnalysisProcessor(samplesdict, dtype=dtype)
    elif proc == 'centralGen':
        import centralGen
        processor_instance = centralGen.AnalysisProcessor(samplesdict, dtype)
    elif proc == 'shellProc':
        import shellProc
        processor_instance = shellProc.AnalysisProcessor(samplesdict, dtype)

    output = runner(flist, treename, processor_instance)

    dt = time.time() - tstart

    if executor == "work_queue":
        print(f'Processed {nevts_total} events in {dt:.2f} seconds ({nevts_total/dt:.2f} evts/sec).')
    elif executor == "futures":
        print(f'Processing time: {dt:.2f} with {nworkers} ({dt*nworkers:.2f} cpu overall)')

    if not os.path.isdir(outname): os.system("mkdir -p %s"%outname)
    out_train_features  = os.path.join(outname,"train_features.p")
#    out_train_fit_coefs = os.path.join(outname,"train_fit_coefs.p")
#    out_test_features  = os.path.join(outname,"test_features.p")
#    out_test_fit_coefs = os.path.join(outname,"test_fit_coefs.p")
    print(f"\nSaving output in {outname}...")
    save(output['train_features'].get(), out_features)
#    save(output['train_fit_coefs'].get(), out_fit_coefs)
#    save(output['test_features'].get(), out_features)
#    save(output['test_fit_coefs'].get(), out_fit_coefs)

    print("Done!")
    
if __name__ == '__main__':
    profile = False
    if profile:
        cProfile.run('main()', filename='profile.txt')
    else:
        main()