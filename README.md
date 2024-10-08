# ttbarEFT
Top quard pair EFT analysis using the Coffea framework and discriminators build using [EFTmva](https://github.com/emcgrady/EFTmva)

## Getting started

If conda or mamba is not already available, download and install your preferred package manager. These instructions will use `micromamba`. 

```
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

The ttbarEFT directory is set up to be installed as a python package. This is accomplished with the following series of commands.

```
git clone https://github.com/TopEFT/ttbarEFT.git
cd ttbarEFT
unset PYTHONPATH
micromamba env create -f environment.yml
micromamba activate ttbarEFT
pip install -e .
```

The `-e` option installs the project in editable mode (i.e. setuptools "develop mode"). If you wish to uninstall the package, you can do so by running `pip uninstall ttbarEFT`. The `topcoffea` package upon which this analysis also depends is not yet available on `PyPI`, so we need to clone the `topcoffea` repo and install it ourselves.

```
cd /your/favorite/directory
git clone https://github.com/TopEFT/topcoffea.git
cd topcoffea
pip install -e .
```
Now all of the dependencies have been installed and the ttbarEFT repository is ready to be used. The next time you want to use it, all you have to do is to activate the environment via 

```
unset PYTHONPATH
micromamba activate ttbarEFT
```

