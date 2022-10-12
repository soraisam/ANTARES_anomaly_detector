import numpy as np
from scipy import stats
import json
import pandas as pd
import os
import argparse
import sys

zero_dm = 1.0E-30
zero_log_dt = -7.0

def abi_method(PB, files, directory, cc):
    """This function parses the dm-dt json files."""
    out={}

    for i in range(len(PB)):
        out[PB[i]]={'x':[], 'y':[], 'ids':[]}
        
    for item in files:
        a = item.replace(".json", "")

        path = directory + item
        with open(path) as f:
            b = json.load(f)        
    
        for key, val in b.items():
            if key not in out.keys():
                continue
            tempo = len(val['dmag']) * [a]
            out[key]['x'] = out[key]['x'] + val['log_dt']
            out[key]['y'] = out[key]['y'] + val['dmag']
            out[key]['ids'] = out[key]['ids'] + tempo
        
    return out



def construct_pdf(X, Y, non_empty_edges_1d, PB):
    """Writes the json file that strips the normalized dmags per log dt bin."""
    if (os.path.isfile(f"./{PB}_logdt_dm.json")):
        print (f"{PB}_logdt_dm.json already exists, so not overwriting. If wishing to regenerate, delete it and try again.")                
    else:
        bin_means, bin_edges, binnumber = stats.binned_statistic(X, Y, statistic='count', bins=non_empty_edges_1d)        
        Big_to_write = {}
        Big_to_write['log_dt_bin_edges'] = list(bin_edges)
        
        for k in range(1, len(bin_edges)):
            values = Y[binnumber==k]
            dm_hist, dm_bin_edges = func_1d_hist(values)
            Big_to_write[str(k)] = {'dm_hist':list(dm_hist), 'dm_bin_edges':list(dm_bin_edges)}
            
        with open(f"{PB}_logdt_dm.json", "w") as f:
            json.dump(Big_to_write, f, sort_keys=True, indent=4)
            print (f"Wrote the P(dm|logdt) file for {PB}")



def get_non_empty_bins(X):
    hist_1d, bin_edges_1d = np.histogram(X, bins='auto') 
    # Here join the empty bins (right merging the bins) such that there is at least 100 points per log-dt bin. 

    non_empty_edges_1d=[bin_edges_1d[0]]
    non_empty_hist_1d=[]
    
    per_bin_min=100

    print (f"Number of bins before grouping empty ones is {len(hist_1d)}")
    i=0
    while (i<len(hist_1d)):
        if (hist_1d[i]>=per_bin_min):
            non_empty_edges_1d.append(bin_edges_1d[i+1])
            non_empty_hist_1d.append(hist_1d[i])
            i = i+1
        else:
            Sum = 0
            for j in range(i,len(hist_1d)):
                Sum = Sum + hist_1d[j]
                if Sum>=per_bin_min or j==(len(hist_1d)-1):
                    non_empty_edges_1d.append(bin_edges_1d[j+1])
                    non_empty_hist_1d.append(Sum)
                    i=j+1
                    break

    print (f"Number of bins after grouping empty ones is {len(non_empty_hist_1d)}")
    return non_empty_edges_1d




def func_1d_hist(dm):
    hist, edges = np.histogram(dm, bins='auto', density=True) 
    mask = hist==0.0
    hist[mask] = zero_dm
    return hist, edges



def main(data_cat_path='./new_jsons/', passband=None, cc=None):
    """
    Parameters
    ----------
    data_cat_path : str
        The relative path of the folder with the dm-dt json files of the reference light curves.
    passband : list
        The unique passbands of the reference light curves.  
    cc : dict
        The cross-passbands to be considered for computing the color-based P(dm|dt) distributions.
        Example for ZTF is '{1:g-R, 2:R-g}'
    """        
    filelists = os.listdir(data_cat_path)
    print (f"{len(filelists)} source(s) in the reference sample")
      
    grand_PB = []
    if passband is not None:
        grand_PB = grand_PB + passband
    if cc is not None:
        for K, V in cc.items():
            #grand_PB.append(V[0]+"-"+V[1])
            grand_PB.append(V)

    pdf_models_new = {}
   
    dmdt_data = abi_method(grand_PB, filelists, data_cat_path, cc)
    
    for i in range(len(grand_PB)):
        print ("\nwriting pdf file", grand_PB[i])
        X = np.array(dmdt_data[grand_PB[i]]['x'])
        Y = np.array(dmdt_data[grand_PB[i]]['y'])
        I = np.array(dmdt_data[grand_PB[i]]['ids'])
        non_empty_edges_1d = get_non_empty_bins(X)
        construct_pdf(X, Y, non_empty_edges_1d, grand_PB[i])        
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument('--data_cat', type=str, help='Path of the folder with the dm-dt json files of the reference light curves', default='./new_jsons/')
    parser.add_argument('--passband', type=str, help='Passbands (e.g., R, g for ZTF)', default=None)
    parser.add_argument('--color_code', type=str, help="Example for ZTF is '{\"1\":\"g-R\", \"2\":\"R-g\"}'. The single and double quotes are required.", default=None)

    args = parser.parse_args()
    if args.passband is None and args.color_code is None:
        sys.exit("Please input one of passband or color_code to run this code. Doing nothing exiting.")
        
    if args.passband is not None:
        passband = list(args.passband.replace(' ','').replace(',',''))
    else:
        passband = args.passband
    if args.color_code is not None:
        color_code = json.loads(args.color_code)
    else:
        color_code = args.color_code
    main(args.data_cat, passband, color_code)
    
