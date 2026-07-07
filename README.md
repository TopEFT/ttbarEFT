# ttbarEFT
Top quard pair production EFT analysis using the Coffea framework and discriminators build using [EFTmva](https://github.com/emcgrady/EFTmva)

## Setup (for coffea v2025.7.3)
If conda or mamba is not already available, download and install your preferred package manager. 

Micromamba does not work to install coffea v2025.7.3. Instead use conda or mamba.
The packaging environment step is ***much*** faster if mamba is used to install the environment instead of conda. 

The ttbarEFT directory is set up to be installed as a python package. This is accomplished with the following series of commands.

```bash
git clone https://github.com/TopEFT/ttbarEFT.git
cd ttbarEFT
unset PYTHONPATH
mamba env create -f environment.yml
mamba activate ttbarEFT
pip install -e .
```

### topcoffea installation
The `-e` option installs the project in editable mode (i.e. setuptools "develop mode"). If you wish to uninstall the package, you can do so by running `pip uninstall ttbarEFT`. The `topcoffea` package upon which this analysis also depends is not yet available on `PyPI`, so we need to clone the `topcoffea` repo and install it ourselves.
```bash
cd 
git clone https://github.com/TopEFT/topcoffea.git
cd topcoffea
pip install -e .
```

### vineReduce installation
Now setup the `vineReduce`  repository, originally named `dynamic_data_reduction`. 
Referencing the instructions from the [vineReduce repo](https://github.com/cooperative-computing-lab/dynamic_data_reduction#installing-from-source), install this repo inside your ttbarEFT environment
```bash
cd 
git clone https://github.com/cooperative-computing-lab/dynamic_data_reduction.git
cd dynamic_data_reduction
pip install -e .
```

### ttbarEFT  installation 
Now all of the dependencies have been installed and the ttbarEFT repository is ready to be used. The next time you want to use it, all you have to do is to activate the environment via 
```bash
unset PYTHONPATH
unset PERL5LIB
mamba activate ttbarEFT
```

# Submit Workers 
### Package Environment
Currently, there is an issue with taskvine such that using `run_analysis.py` where the remote environment file is passed using the `environment_file` option fails with `ModuleNotFoundError` for `cloudpickle`. However, it works to provide the packaged environment tarball to the taskvine worker directly. 

To make the packaged env, use the script in this analysis directory to run `remote_environment.py`. 

Note: This script calls the version of `remote_enironment.py` that is in this repository, not the `topcoffea` version. Since September 2025, the `topcoffea` version pins old versions of coffea, setuptools, and ndcctools which are not compatible with coffea2025. 

Run the script with: 
```bash
unset PYTHONPATH
unset PERL5LIB
mamba activate ttbarEFT

cd ttbarEFT/analysis 
python make_packaged_env.py
```

Take note of the tarball path. It should be something like `ttbarEFT/analysis/topeft-envs/env_spec_*_HEAD.tar.gz`. This will be used in the commands below to submit workers. 

### Submit workers on Tier3 and/or workers to condorfe which don't use xrootd 
Then, you can submit workers to the Tier3 or condorfe using this command: 
```bash
unset PYTHONPATH
unset PERL5LIB
mamba activate ttbarEFT

vine_submit_workers -T condor -M ${USER}-taskvine-coffea --python-env path/to/tarball/env -t 900 --cores Ncores --memory Nmemory --disk disk Nworkers
```

### Submit workers to condorfe which use xrootd 
To submit workers to condorfe that access files on xrootd, first active your voms proxy, and then use this command to include the `with_oasis_certs` script in order to set the permissions for xrootd correctly. 
In this case, `with_oasis_certs` is in the same directory that this command  is being run from. Update accordingly.  
```bash
unset PYTHONPATH
unset PERL5LIB
mamba activate ttbarEFT

vine_submit_workers --extra-file with_oasis_certs --wrapper ./with_oasis_certs -T condor -M ${USER}-ddr-coffea --python-env path/to/tarball/env --cores Ncores --memory Nmemory --disk disk Nworkers
```

# Making and Updating Sample Json Files 
1. Edit `ttbarEFT/scripts/make_jsons.py` to add the samples you want to produce jsons for.
2. Run this script `ttbarEFT/scripts/python make_jsons.py` (this calls `ttbarEFT/scripts/make_sample_json.py` in a loop).
3. Copy the sample jsons to the appropriate directory in `input_samples/sample_jsons/`.
4. Make a config (yaml) that points to the new jsons and put it in `input_samples/cfgs/`.
5. Run `update_sow_jsons.py` using the processor `sow_processor.py` over the new config. This script gets the up/down variations for ISR/FSR/renorm/fact and also corrects the sum of weights for EFT samples.
```bash
unset PYTHONPATH
unset PERL5LIB
mamba activate ttbarEFT

cd ttbarEFT/analysis/
python update_sow_jsons.py -p sow_processor.py -x ddr --relpath '../input_samples/sample_jsons/background_samples/central_UL/' -o <output file name> <path to config>
```
6. If you want to create a yaml for a skimmed sample, first do step #5. The sow variations have to be calculated on the unskimmed sample so those need to be available in the original json file. Then, run the script `ttbarEFT/scripts/update_jsons.py` (a wrapper for `ttbarEFT/scripts/update_sow_jsons.py`) for many jsons at once, or run `update_sow_jsons.py` for individual sample jsons. This script reads in the existing json and copies over quantities that don't change (e.g. xsec, year, dataset, histAxisName, WC, sow + variations) and updates the file names to point to files in the directory you provide. The result is a new sample json for a skimmed sample where just the file list is updated. 

# Repository Contents
Example to run the main processor over all Run2 samples: 
```bash
python run_processor_vineReduce.py -p analysis_processor.py --syst-list "all" --hist-list mllbb ptll l0pt -x ddr --doSR --doPDF -o SR_ALL_2j3j --chunksize 250000 ../input_samples/cfgs/Run2_all_MCandDATA.yaml
```


`ttbarEFT/analysis`
- `analysis_process.py`: Main analysis processor! Includes all selections, corrections, etc. The SR channels dictionary can be changed in the `__init__` in the lines like `self._channels = cat_dict['SR_CHANNELS_mllbb'][lep_cat]`. This way, the SR dictionary selected can be changed without remaking the packaged env. 
- `make_CR_plots.py`: Main plotting script. This has wrapper functions to make CR plots, SR plots, LOtoNLO comparison plots, etc. 
- `run_sow.py`& `sow_processor.py` : processor and corresponding run script to get sow with variations. 
- `update_sow_jsons.py`: Alternative to `run_sow.py`. This updates the json files with the sow found in `sow_processor.py`. See instructions above. 
- `get_sow.py`: This script reads the pkl produced by `sow_processor.py` and prints out the values
- `make_packaged_env.py`: Script to package up the conda environment. (See instructions above).


- `LOtoNLOuncert.ipynb`: notebook to derive the missing syst uncertainty to cover the gap between SMEFT LO sample + theory uncert and Powheg NLO sample. Produces a json that is copied to `ttbarEFT/params/ttLOuncert.json`
- `Neff_processor.py`: processor that makes histograms with the sumw and sumw2 - used to calculate the effective number of events in the original SMEFTsim samples to figure out how many events were needed in the recreated sample

- Notebooks used for quick checks: 
    - `CheckProcessYields.ipynb`: notebook to return the MC predicted yield from a pkl file - used to extract predicted yields to compare against what the datacard maker produced 
    - `MET_checks.ipynb`: used to confirm that the MET unclustered energy variation was working
    - `checkEFTeffects.ipynb`: notebook to check the strength of EFT effects to determine which bins were safe to look at before fully unblinding
    - `unblinding.ipynb`: notebook to check the first few bins unblinded to try and catch any major mistakes (based on `checkEFTeffects` notebook)  
- Notebooks used for development of features/corrections/etc. All of these should not be needed in the future since the final working versions of these features are all implemented in the main analysis scripts. 
    - `PDFUncertainties.ipynb`: used for initial production/tests of the PDF variation
    - `dev_JEC.ipynb`: used while developing the jet energy corrections
    - `dev_muontrig_uncert.ipynb`: fixing and testing the error propagation of muon trigger uncertainty
- Old version of processors and plotting scripts that can probably be deleted 
    - `make_CR_plots_ddr.py`
    - `make_CR_plots_origddr.py`
    - `optimize_analysis_processor.py`
    - `optimize_processor_noJER.py`
    - `optimize_run_processor_vineReduce.py`
    - `refactored_run_processor_with_ddr.py`
    - `run_processor_vineReduce.py`

`ttbarEFT/LHEmtt`
- Contains the scripts to find the mtt stitching factor 

`ttbarEFT/LOtoNLO`
- scripts to find the LOtoNLO reweighting
- this version of `make_CR_plots.py` is setup to make comparisons of Powheg vs uncorrected SMEFT vs corrected SMEFT
- `ttbarLOtoNLO_reweighting.ipynb` notebook used to develop the reweighting factor
- `make_reweighting_processor.py` & `run_reweight_processor.py` : processor and run script to create the pkl files to derive the LOtoNLO reweighting. Example usage below
- `apply_reweighting_processor.py` & `run_processor_vineReduce.py`: processor and run script to apply the reweighting to make comparison plots between SMEFT and Powheg. These correctly normalize the sum of weights. 
examples of run commands: 
```bash
python run_reweight_processor.py -x ddr -p make_reweighting_processor.py --chunksize 500000 -o SMEFTsim_avgpt_smallerbins_var2  ../input_samples/cfgs/Run2_TT01j2l_EFT.yaml

python run_reweight_processor.py -x ddr -p make_reweighting_processor.py --chunksize 500000 -o Powheg_avgpt_smallerbins_var2 ../input_samples/cfgs/Run2_PowhegTTTo2L2Nu.yaml
```

```bash
python run_processor_vineReduce.py -p apply_reweighting_processor.py --syst-list "all" --hist-list mllbb l0eta l0phi l0pt j0pt j0eta j0phi ptllbb ptll lldphi -x ddr -o SMEFTsim_toppt_rwgt --doSR ../input_samples/cfgs/2017_TT01j2l_EFT.yaml
```

`ttbarEFT/btagMCeff`
- Scripts to derive the btagging MC eff 2d histograms. 
- `btagMCeff_processor.py` & `run_btagMCeff_processor.py`: processor and run script to make the btagMCeff histograms. The central Powheg sample is used after it was confirmed that the btagMCeff are basically the same when derived using Powheg, SMEFTsim. The produced pkl for each era is copied into `ttbarEFT/data/btagSF/UL` and then used by the analysis_processor. 
- `make_btagMCeff_hists.py`: plots the btag histograms

example usage (repeat for other eras)
```bash
python run_btagMCeff_processor.py -p btagMCeff_processor.py -o btagMCeff_addJEC_2016APV -x ddr ../input_samples/cfgs/btag_eff.yaml

python make_btagMCeff_hists.py --file btagMCeff_2016APV.pkl.gz --outdir btagMCeff_addJEC_2016APV/
```

`ttbarEFT/mc_validation`: 
- A few scripts copied over from the old MC validation I did. But all of those old scripts were not setup to run with coffea2025. These are not fully fixed and would need some updating but the basic idea is there.

`ttbarEFT/datacard_maker`
- makes the datacards for Combine. Basically a wrapper for `ttbarEFT/modules/datacard_tools.py`
example usage
```bash
python make_cards.py "../analysis/SR_ALL_2j3j_260519.pkl.gz" --out-dir "260520_FINAL" --do-mc-stat --do-nuisance --unblind
```

# Combine
Everything is based on the [py3Combine branch of the EFTFit repo](https://github.com/TopEFT/EFTFit/tree/py3Combine). I didn't try the "new fancy install script". A summary of the setup instructions that I used is below but this should be identical to that README.
Then checkout ttbarEFT. 
- The main two scripts (`Fitter/scripts/EFTFitter.py` and `Fitter/scripts/EFTPlotter.py`) have a lot of hard-coding for the WCs used in the multilepton analysis. 
- The versions in the `ttbarEFT` branch are updated for this analysis. 
- Also I updated the script to get the whole list of nuisance parameters from the workspace instead of having to hard-code that inside `EFTFitter.py`
- Unfortunately the scripts are really hard-coded to run inside `Fitter/test` for submitting to condor so if you want to run multiple times, you have to rename or remove the existing `test` dir and then rerun

## Setup OUTSIDE of any conda env: 
```bash
export SCRAM_ARCH=el9_amd64_gcc12
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv

git clone git@github.com:cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit/
git fetch origin
git checkout v10.0.2
cd -
scramv1 b clean; scramv1 b

cd $CMSSW_BASE/src/
git clone https://github.com/TopEFT/EFTFit.git EFTFit
git checkout ttbarEFT
scramv1 b clean; scramv1 b

git clone git@github.com:cms-analysis/CombineHarvester.git
cd CombineHarvester
git checkout 128e41eb
scramv1 b clean; scramv1 b
```

## Running Combine
This is assuming that you already ran the `make_cards.py` script to create datacards from the pkl. Change paths as needed. 

### Make workspace: 
```bash
cd ~/CMSSW_14_1_0_pre4/src/

unset PYTHONPATH
unset PERL5LIB
cmsenv

cd EFTFit/Fitter/
mkdir test/
cd test/

cp ~/ttbarEFT-coffea2025/datacard_maker/260520_FINAL/*.root . 
cp ~/ttbarEFT-coffea2025/datacard_maker/260520_FINAL/*.txt .
cp ~/ttbarEFT-coffea2025/datacard_maker/260520_FINAL/scalings.json .  
cp ~/ttbarEFT-coffea2025/datacard_maker/260520_FINAL/selectedWCs.txt .

combineCards.py TESTING-ttbar*.txt > combinedcard.txt
time source ../scripts/text2workspace.sh
```

### Run 1d scans:
To unblind, remove `other=['-t', '-1']`

Change the basename to what you want the titles to be. The exact basename used here then has to be used for the plotting. 
```bash
python3 -i ../scripts/EFTFitter.py

profiled: 
fitter.batch1DScanEFT(basename='260520_FINAL_asimov-profiled', batch='condor', workspace='workspace.root', other=['-t', '-1'])

frozen: 
fitter.batch1DScanEFT(basename='260520_FINAL_asimov-frozen', batch='condor', workspace='workspace.root', freeze=True, other=['-t', '-1'])

fitter.batchRetrieve1DScansEFT(basename='260520_FINAL_asimov-profiled', batch='condor')
fitter.batchRetrieve1DScansEFT(basename='260520_FINAL_asimov-frozen', batch='condor')

python3 -i ../scripts/EFTPlotter.py
plotter.BestScanPlot(basename_float_lst='260520_FINAL_asimov-profiled', basename_freeze_lst='260520_FINAL_asimov-frozen', filename='260520_FINAL_asimov', titles=['Others Profiled', 'Others Frozen'], printFOM=True)

plotter.BatchOverlayLLPlot1DEFT(basename1_lst=['260520_FINAL_asimov-profiled'], basename2_lst=['260520_FINAL_asimov-frozen'], wcs=[], log=False, final=False, titles=['Others profiled', 'Others fixed to SM'])
```

### Summary Plots 
```bash
python3 -i ../scripts/EFTPlotter.py
plotter.BestScanPlot(basename_float_lst='260520_FINAL_asimov-profiled', basename_freeze_lst='260520_FINAL_asimov-frozen', filename='260520_FINAL_asimov', titles=['Others Profiled', 'Others Frozen'], printFOM=True)

plotter.BatchOverlayLLPlot1DEFT(basename1_lst=['260520_FINAL_asimov-profiled'], basename2_lst=['260520_FINAL_asimov-frozen'], wcs=[], log=False, final=False, titles=['Others profiled', 'Others fixed to SM'])
```

`BestScanPlot` makes the summary plot and also prints out the 1$\sigma$ and 2$\sigma$ CIs and best fit values. 

`BatchOverlayLLPlot1DEFT` make the 1d scan plots with the profiled and frozen both included. 

### Run impact plots
blinded: 
```bash
python3 -i ../scripts/EFTFitter.py
fitter.ImpactInitialFit(workspace='workspace.root', wcs=[])
fitter.ImpactNuisance(workspace='workspace.root', wcs=[])
fitter.ImpactCollect(workspace='workspace.root', wcs=[])
```

unblinded: 
```bash
python3 -i ../scripts/EFTFitter.py
fitter.ImpactInitialFit(workspace='workspace.root', unblind=True, wcs=[])
fitter.ImpactNuisance(workspace='workspace.root', unblind=True, wcs=[])
fitter.ImpactCollect(workspace='workspace.root', unblind=True, wcs=[])
```

### Optionally use more robust fitting settings for ctG: 
```bash
python3 -i ../scripts/EFTFitter.py

fitter.batch1DScanEFT_ctG(basename='.testRandProf_prof_ctG', batch='condor', workspace='workspace.root', RandProf = 25, scan_wcs=['ctGRe', 'ctGIm'])
fitter.batch1DScanEFT_ctG(basename='.testRandProf_froz_ctG', batch='condor', workspace='workspace.root', RandProf = 25, scan_wcs=['ctGRe', 'ctGIm'], freeze=True)

fitter.batchRetrieve1DScansEFT(basename='.testRandProf_prof_ctG', batch='condor')
fitter.batchRetrieve1DScansEFT(basename='.testRandProf_froz_ctG', batch='condor', scan_wcs=['ctGRe', 'ctGIm'])
```

```bash
python3 -i ../scripts/EFTPlotter.py
plotter.BatchLLPlot1DEFT(basename_lst=['.testRandProf_prof_ctG'], wcs=['ctGRe', 'ctGIm'], log=False)
```

# Misc. Notes
- The version of the run script and analysis processor for usage with Coffea Executors (pre-refactoring) has been moved to `coffeaExector` so it's available if testing is needed.

# Previous Analysis
- Copy of TOP-17-020: https://arxiv.org/pdf/1903.11144
- Reza's NanoAOD analysis script: https://github.com/rgoldouz/EFTNanoAnalysis/blob/master/NanoAnalysis/src/MyAnalysis.cc
- Reza's FCNC analysis repo: https://github.com/rgoldouz/FCNC/tree/main/NanoAnalysis
