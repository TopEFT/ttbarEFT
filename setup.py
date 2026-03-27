import setuptools

setuptools.setup(
    name='ttbarEFT', 
    version='0.0.0',
    description = 'Analysis code for ttbar EFT analysis',
    packages=setuptools.find_packages(),
    package_data={
        "ttbarEFT": [
            "params/*.json",
            "params/*.yaml",
            "data/POG/*/*.json.gz",
            "data/POG/JME/*/*.json.gz",
            "data/POG/EGM/*/*.json.gz",
            "data/POG/MUO/*/*.json.gz",
            "data/POG/MUO/*/*.json",
            "data/POG/LUM/*/*.json.gz",
            "data/POG/BTV/*/*.json.gz",
            "data/btagSF/UL/*",
            "data/leptonSF/elec/*.root",
            "modules/*.yaml",
        ]
    }
)
