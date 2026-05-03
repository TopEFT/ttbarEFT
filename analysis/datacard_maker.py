import numpy as np

import cloudpickle
import argparse
import uproot
import gzip
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--discriminator', 
                       default = '.', 
                       help = 'path to the directory containing .gz file containing the discriminator.')
    parser.add_argument('--root', '-r', default='dnn_disc', help='name of the output template root file')

    args = parser.parse_args()
    
    out_root = args.root
    prefix   = args.discriminator
    
    file   = f'{prefix}/histeft.pkl.gz'
    
    with gzip.open(file,"rb") as f:
        histo = cloudpickle.load(f)

    histo = histo['genDNNDisc']
    with uproot.recreate(f"{prefix}/dnn_disc.root") as file:
        file["data_obs"] = histo[{'cat': 'genDNNDisc'}].as_hist({'ctq8': 0})
        file["ttbar"] = histo[{'cat': 'genDNNDisc'}].as_hist({'ctq8': 0})
        
    scaling = histo[{'cat': 'genDNNDisc'}].make_scaling()
    scaling[-2] += scaling[-1]
    scaling[1] += scaling[0]
    scaling = scaling[1:-1]
    mask = scaling[:,0] != 0
    scaling[mask] /= np.expand_dims(scaling[mask,0], 1)
    
    out_json = [
        { 
            "channel"    : "bin_1",
            "process"    : "ttbar",
            "parameters" : [
                "cSM[1]",
                "ctq8[0,-10,10]"
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
    num     = 0
    
    bin_str = f"bin_1"
    col_width = 5
    sgnl_count = histo.as_hist({'ctq8': 0}).sum()
    
    with open(outf_card_name,"w") as f:
        f.write(f"imax 1\n")
        f.write(f"jmax 0\n")
        f.write(f"kmax 0\n")
        f.write(line_break)
        f.write(f"shapes *        * dnn_disc.root $PROCESS\n")
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

        row = [f"{'process':<{left_width}}"]
        row.append(f"{num:>{col_width}}")
        row = " ".join(row) + "\n"
        f.write(row)
        
        row = [f"{'rate':<{left_width}}"]
        row.append(f"{sgnl_count:>{col_width}.2f}")
        row = " ".join(row) + "\n"
        f.write(row)

        f.write('* autoMCStats 10')

if __name__ == '__main__':
    main()