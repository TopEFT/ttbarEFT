import topcoffea.modules.remote_environment as remote_environment

if __name__ == '__main__':

	remote_environment.get_environment(
                extra_pip_local = {"ttbarEFT": ["ttbarEFT", "setup.py"]},
                # extra_pip=['mt2'],
                # extra_conda=["pytorch=2.3.1", "numpy=1.23.5"]
			)