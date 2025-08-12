import setuptools

setuptools.setup(
    name='ttbarEFT', 
    version='0.0.0',
    description = 'Analysis code for ttbar EFT analysis',
    packages=setuptools.find_packages(),
    package_data={
        "ttbarEFT": [
            "params/*.json",
            "data/POG/JME/*/*.json.gz"
        ]
    }
)
