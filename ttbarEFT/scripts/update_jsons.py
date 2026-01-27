from update_sample_json_files import update_json
import os

def update_multiple_jsons(ref_dir, path_dir, json_dict):
	
	for refname in json_dict.keys():
		ref = os.path.join(ref_dir, refname)
		path = os.path.join(path_dir, json_dict[refname])
		outname = f"{json_dict[refname]}_ptSkim"

		if os.path.exists(f"{outname}.json"): 
			print(f"{outname} file already exists. Skipping...")
			continue
		else: 
			update_json(path=path, refpath=ref, outname=outname)

# user inputs
data_ref_dir = "/users/hnelson2/ttbarEFT-coffea2025/input_samples/sample_jsons/data_samples/2017/"
data_path_dir = "/cms/cephfs/data/store/user/hnelson2/skims/data/FullRun2/ptSkim/"

# 'reference json': 'directory with new files'
data_2017_dict = {
	'DoubleEG_B-UL2017.json': 'DoubleEG_B_UL2017',
	'DoubleEG_C-UL2017.json': 'DoubleEG_C_UL2017',
	'DoubleEG_D-UL2017.json': 'DoubleEG_D_UL2017',
	'DoubleEG_E-UL2017.json': 'DoubleEG_E_UL2017',
	'DoubleEG_F-UL2017.json': 'DoubleEG_F_UL2017',
	'SingleMuon_B-UL2017.json': 'SingleMuon_B_UL2017',
	'SingleMuon_C-UL2017.json': 'SingleMuon_C_UL2017',
	'SingleMuon_D-UL2017.json': 'SingleMuon_D_UL2017',
	'SingleMuon_E-UL2017.json': 'SingleMuon_E_UL2017',
	'SingleMuon_F-UL2017.json': 'SingleMuon_F_UL2017',
}

bkgd_ref_dir = "/users/hnelson2/ttbarEFT-coffea2025/input_samples/sample_jsons/background_samples/central_UL/"

DY_path_dir = "/cms/cephfs/data/store/user/hnelson2/skims/mc/ptSkim/DY/"
DY_dict = {
	'UL16APV_DY10to50.json'	: 'UL16APV_DY10to50',
	'UL16APV_DY50.json'		: 'UL16APV_DY50',
	'UL16_DY10to50.json'	: 'UL16_DY10to50',
	'UL16_DY50.json'		: 'UL16_DY50',
	'UL17_DY10to50.json'	: 'UL17_DY10to50',
	'UL17_DY50.json'		: 'UL17_DY50',
	'UL18_DY10to50.json'	: 'UL18_DY10to50',
	'UL18_DY50.json'		: 'UL18_DY50',
}

bkgd_path_dir = "/cms/cephfs/data/store/user/hnelson2/skims/mc/ptSkim/background/"
bkgd_2017_dict = {
	'UL17_TTGJets.json'			: 'UL17_TTGJets',
	'UL17_TTWJetsToLNu.json'	: 'UL17_TTWJetsToLNu',
	'UL17_TTZToLLNuNu_M_10.json': 'UL17_TTZToLLNuNu_M_10',
	'UL17_WJetsToLNu.json'		: 'UL17_WJetsToLNu',
	'UL17_WWTo2L2Nu.json'		: 'UL17_WWTo2L2Nu',
	'UL17_WWW_4F.json'			: 'UL17_WWW_4F',
	'UL17_WWZ_4F.json'			: 'UL17_WWZ_4F',
	'UL17_WZTo3LNu.json'		: 'UL17_WZTo3LNu',
	'UL17_WZZ.json'				: 'UL17_WZZ',
	'UL17_ZZTo4L.json'			: 'UL17_ZZTo4L',
	'UL17_ZZZ.json'				: 'UL17_ZZZ',
}

update_multiple_jsons(data_ref_dir, data_path_dir, data_2017_dict)
update_multiple_jsons(bkgd_ref_dir, DY_path_dir, DY_dict)
update_multiple_jsons(bkgd_ref_dir, bkgd_path_dir, bkgd_2017_dict)