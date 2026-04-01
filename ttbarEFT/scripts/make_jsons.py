from make_sample_json import make_sample_json
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