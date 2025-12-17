# ttbarEFT
Top quard pair EFT analysis using the Coffea framework and discriminators build using [EFTmva](https://github.com/emcgrady/EFTmva)

## Setup (for coffea v2025.7.3)

If conda or mamba is not already available, download and install your preferred package manager. 
Micromamba does not work to install coffea v2025.7.3. Instead use conda (tested) or mamba (untested)

The ttbarEFT directory is set up to be installed as a python package. This is accomplished with the following series of commands.

```
git clone https://github.com/TopEFT/ttbarEFT.git
cd ttbarEFT
unset PYTHONPATH
conda env create -f environment.yml
conda activate ttbarEFT
pip install -e .
```

The `-e` option installs the project in editable mode (i.e. setuptools "develop mode"). If you wish to uninstall the package, you can do so by running `pip uninstall ttbarEFT`. The `topcoffea` package upon which this analysis also depends is not yet available on `PyPI`, so we need to clone the `topcoffea` repo and install it ourselves.

```
cd 
git clone https://github.com/TopEFT/topcoffea.git
cd topcoffea
pip install -e .
```
Now all of the dependencies have been installed and the ttbarEFT repository is ready to be used. The next time you want to use it, all you have to do is to activate the environment via 

```
unset PYTHONPATH
unset PERL5LIB
conda activate ttbarEFT
```

## Submit Workers 
Currently, there is an issue with taskvine such that using `run_analysis.py` where the remote environment file is passed using the `environment_file` option failes with `ModuleNotFoundError` for `cloudpickle`. However, providing the packaged environment tarball to the taskvine worker succeedes. 


**First**, 
Update your copy of `topcoffea/topcoffea/modules/remote_environment.py`. This file currently (September 11, 2025) pins old versions of coffea , setuptools, and ndcctools. These need to be udpated to work with coffea v2025.7.3. Modify lines 25-38 to contain the following: 

```
default_modules = {
    "conda": {
        "channels": ["conda-forge"],
        "packages": [
            f"python={py_version}",
            "awkward=2.8.7",
            "coffea=2025.7.3",
            "numpy",
            "ndcctools",
            "pip",
            "conda",
            "conda-pack",
            "dill",
            "xrootd",
            "setuptools==80.9.0",
        ],
    },
    "pip": ["topcoffea"],
}
```


**Second**, create the packaged conda environment using 
```
unset PYTHONPATH
unset PERL5LIB
conda activate ttbarEFT

cd ttbarEFT/analysis 
python make_packaged_env.py
```


**Third**, make sure that in your version of `run_processor.py`, `environment_file` in `executor_args` is **commented out**. 
Take note of the tarball path. It should be something like `ttbarEFT/analysis/topeft-envs/env_spec_*_HEAD.tar.gz`.


**Forth**, in a separate terminal instance, submit workers using this command 
```
unset PYTHONPATH
unset PERL5LIB
conda activate ttbarEFT

vine_submit_workers -T condor -M ${USER}-taskvine-coffea --python-env <PathToTarball> -t 900 --cores 4 --memory 16000 --disk 100000 10
```
## Setup (for running with VineReduce) 
dynamic data reduction repo: [repo](https://github.com/cooperative-computing-lab/dynamic_data_reduction)

example usage: 
```bash
python run_processor_vineReduce.py -x ddr -p analysis_processor_ddr.py -o MC_CRs_ddr ../input_samples/cfgs/MC_CR_2017.yaml
```

submit workers: 
```
vine_submit_workers --extra-file with_oasis_certs --wrapper ./with_oasis_certs -T condor -M ${USER}-ddr-coffea --python-env /users/hnelson2/ttbarEFT-coffea2025/analysis/topeft-envs/env_spec_251e393e_edit_HEAD.tar.gz --cores 12 --memory 16000 --disk 300000 20
```

# Notes
The version of the run script and analysis processor for usage with Coffea Executors (pre-refactoring) has been moved to `coffeaExector` so it's available if testing is needed.


# Previous Analysis

- Copy of TOP-17-020: https://arxiv.org/pdf/1903.11144
- Reza's NanoAOD analysis script: https://github.com/rgoldouz/EFTNanoAnalysis/blob/master/NanoAnalysis/src/MyAnalysis.cc
- Reza's FCNC analysis repo: https://github.com/rgoldouz/FCNC/tree/main/NanoAnalysis
