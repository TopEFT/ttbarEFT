import correctionlib as clib

def get_available_corrections(json_path):
    cset = clib.CorrectionSet.from_file(json_path)
    available_corrections = list(cset.keys())

    return available_corrections


# print(get_available_corrections("/users/hnelson2/topcoffea/topcoffea/data/POG/JME/2022_Summer22/jet_jerc.json.gz"))

print(get_available_corrections("/users/hnelson2/ttbarEFT-coffea2025/ttbarEFT/data/POG/JME/2016preVFP_UL/jet_jerc.json.gz"))