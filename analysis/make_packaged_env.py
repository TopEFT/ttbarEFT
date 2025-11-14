# import topcoffea.modules.remote_environment as remote_environment
import ttbarEFT.modules.remote_environment as remote_environment

if __name__ == '__main__':

	remote_environment.get_environment(
                # extra_pip_local = {"ttbarEFT": ['ttbarEFT', 'setup.py'],},
                extra_pip_local = {"ttbarEFT": ["ttbarEFT", "setup.py"],
                                     "dynamic_data_reduction": []},
                # extra_pip=['mt2'],
                extra_conda=["fsspec-xrootd"]
			)
