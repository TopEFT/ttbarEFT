import numpy as np

import cloudpickle
import argparse
import uproot
import gzip
import json

def get_wc_names_cross(wc_names_lst):
    '''
    returns the list of names for WCs and cross terms
    '''
    wcs = {}
    index = 0
    for i in range(len(wc_names_lst)):
        for j in range(i+1):
            wcs[(wc_names_lst[i], wc_names_lst[j])] = index
            index += 1 

    return wcs

def make_scaling(self, wc_list=None):
    '''
    returns scalings list with flow bins
    wc_list: list of WCs if different than the list contained in the passed HistEFT
    ''' 
    
    if wc_list is None:
        wcs = get_wc_names_cross(['sm'] + self.wc_names)
        scaling = self.values(flow=True)[:,1:-1]

    else:
        wcs =  get_wc_names_cross(['sm'] + wc_list)
        old_wcs = get_wc_names_cross(['sm'] + self.wc_names)
        old_scaling = self.values(flow=True)[:,1:-1]
        scaling = np.zeros((old_scaling.shape[0], len(wcs)))
        for key in wcs.keys():
            if key in old_wcs.keys():
                scaling[:,wcs[key]] = old_scaling[:,old_wcs[key]]
    
    scaling = (scaling/np.expand_dims(scaling[:,0], 1))
    
    for key in wcs.keys():
        if key[0] != key[1]:
            scaling[:,wcs[key]] /= 2
            
    return scaling

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--discriminator', 
                       default = '/scratch365/cmcgrad2/post_processing/net_disc', 
                       help = 'path to the directory containing .gz file containing the discriminator.')
    parser.add_argument('--root', '-r', default='dnn_disc', help='name of the output template root file')

    args = parser.parse_args()
    
    out_root = args.root
    prefix   = args.discriminator
    
    file   = f'{prefix}/histeft.pkl.gz'
    
    with gzip.open(file,"rb") as f:
        histo = cloudpickle.load(f)

    histo = histo['disc']
    file = uproot.recreate(f"{prefix}/dnn_disc.root")
    file['disc'] = histo.as_hist({'ctq8': 0})

    scaling = make_scaling(histo[{'cat': 'net_disc'}])
    
    out_json = [
        { 
            "channel"    : "ch1",
            "process"    : "ttbar",
            "parameters" : [
                "ctq8[0,-1.4,1.4]"
            ],
            "scaling": scaling.tolist()
        }
    ]

    with open(f'{prefix}/scalings.json', 'w') as f:
        json.dump(out_json, f)

    outf_card_name = f'{prefix}/ttbar_EFTmva.txt'
    line_break = "##----------------------------------\n"
    left_width = len(line_break) + 2
    process = 'ttbar'
    
    bin_str = f"bin_1"
    col_width = 5
    sgnl_count = histo.as_hist({'ctq8': 0}).sum()
    
    with open(outf_card_name,"w") as f:
        f.write(f"imax 1\n")
        f.write(f"jmax 0\n")
        f.write(f"kmax 0\n")
        f.write(line_break)
        f.write(f"shapes *        * dnn_disc.root \n")
        f.write(line_break)
        f.write(f"bin         {bin_str}\n")
        f.write(f"observation {sgnl_count:.2f}\n")
        f.write(line_break)
        f.write(line_break)
    
        row = [f"{'bin':<{left_width}}"]
        row.append(f"{bin_str:>{col_width}}")
        row = " ".join(row) + "\n"
        f.write(row)
    
        row = [f"{'process':<{left_width}}"]
        row.append(f"{process:>{col_width}}")
        row = " ".join(row) + "\n"
        f.write(row)
        
        row = [f"{'rate':<{left_width}}"]
        row.append(f"{sgnl_count:>{col_width}.2f}")
        row = " ".join(row) + "\n"
        f.write(row)

if __name__ == '__main__':
    main()