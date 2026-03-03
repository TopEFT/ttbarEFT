import topcoffea.modules.remote_environment as remote_environment
# import ttbarEFT.modules.new_remote_environment as remote_environment
import time

if __name__ == '__main__':

    tstart = time.time()
    remote_environment.get_environment(
                # extra_pip_local = {"ttbarEFT": ['ttbarEFT', 'setup.py'],},
                extra_pip_local = {"ttbarEFT": ["ttbarEFT", "setup.py"],
                                    "dynamic_data_reduction": ["src", "pyproject.toml"]},
                # extra_pip=['mt2'],
                extra_conda=["fsspec-xrootd"],
                # quick=True,
            )

    tend = time.time()
    print(f"\n\n Total time: {tend-tstart}")
