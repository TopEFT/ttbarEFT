from topcoffea.modules.update_json import update_json
import json
import os

def read_json_file(filename):
    with open(filename) as f:
        return json.load(f)


def update_jsons(year, base_path, json_list):
	for name in json_list:
	    orig_json_path = os.path.join(base_path, f"{year}_{name}.json")
	    new_json_path = os.path.join(base_path, f"{year}_{name}_ptSkim.json")
	    
	    if os.path.exists(new_json_path): 
		    json_dict = read_json_file(orig_json_path)
		    keys_to_update = [k for k in json_dict.keys() if k.startswith("nSumOfWeights_")]
		    
		    updates = {i:json_dict[i] for i in keys_to_update}    
		    update_json(new_json_path,dry_run=False,verbose=True, **updates)
	    else: 
	    	print(f"Skipping {new_json_path}. This json doesn't exist.")


background_base_path = "/users/hnelson2/ttbarEFT-coffea2025/input_samples/sample_jsons/background_samples/central_UL/"
background_json_list = [
    'DY10to50',
    'DY50',
    'TTGJets',
    'TTWJetsToLNu',
    'TTZToLLNuNu_M_10',
    'TW_antitop_5f_NoFullyHadronicDecays',
    'TW_top_5f_NoFullyHadronicDecays',
    'WJetsToLNu',
    'WWTo2L2Nu',
    'WWW_4F',
    'WWZ_4F',
    'WZTo3LNu',
    'WZZ',
    'ZZTo4L',
    'ZZZ',
    'TTTo2L2Nu',
]


signal_base_path = "/users/hnelson2/ttbarEFT-coffea2025/input_samples/sample_jsons/signal_samples/"
signal_json_list = [
	'TT01j2l_mtt0to700',
	'TT01j2l_mtt700to900',
	'TT01j2l_mtt900toInf'
]

for y in ['UL16APV', 'UL16', 'UL17', 'UL18']:
	update_jsons(year=y, base_path=background_base_path, json_list=background_json_list)
	# if the signal samples are ever skimmed, this is ready to go
	# # update_jsons(year=y, base_path=signal_base_path, json_list=signal_json_list)

