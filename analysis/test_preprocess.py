import argparse
import json
import time
import cloudpickle
import gzip
import os

import numpy as np
from coffea import dataset_tools

def read_file(filename):
    with open(filename) as f:
        return json.load(f)

if __name__ == '__main__':
	jsonFiles = '/users/hnelson2/ttbarEFT-coffea2025/input_samples/sample_jsons/tests/test.json'

	sample = read_file(jsonFiles)

	print(f"json file contents: {sample}")
	sample_dict = {'test': sample}

	files_dict = {}
	for f in sample['files']: 
		files_dict[f] = {'object_path': 'Events'}

	preprocess_dict = {
		'test': {
			'files': files_dict,
		}
	}

	print(f"\n\n data_dict: \n {preprocess_dict}")

	dataset_runnable = dataset_tools.preprocess(preprocess_dict)[0]
	print(f"\n dataset_runnable: {dataset_runnable}")

	print(f"\n\n USING TWO VARIABLES FROM PREPROCESS \n\n")

	dataset_runnable, dataset_all = dataset_tools.preprocess(preprocess_dict)
	print(f"\n dataset_runnable: {dataset_runnable}")
	print(f"\n dataset_all: {dataset_all}")