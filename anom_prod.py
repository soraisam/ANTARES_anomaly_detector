from antares_client import StreamingClient
import numpy as np
import json
import dmdt_model as dmdt_model
import time 
import argparse
import sys

def get_score_normalized(test_class, pdf_models, expected_D):
    """Computes the expectation-normalized consistency score of a given light curve.
    
    Parameters
    ----------
    test_class : dict of {str : dict}
        The dm-dt values for a given source.    
    pdf_models : dict of {str : dict}
        The conditional probability distributions P(dm|dt).
    expected_D : dict of {str : dict}
        The expectation values of P(dm|dt). 
    """
    
    #zero_dm_log = -2000.0 
    score = {}
    for K, V in test_class.items():
        Sum = 0.0
        for key, val in pdf_models.items():                
            X = np.array(V[key]["x"])
            Y = np.array(V[key]["y"])
            
            mask = Y==-9999    #missing information as opposed to out of range
            X = X[~mask]
            Y = Y[~mask]
            
            if (len(Y)<1):
                continue
                            
            idx_i = np.digitize(X, val['log_dt_bin_edges'])
            
            if 0 in idx_i or (len(val['log_dt_bin_edges']) in idx_i): 
                #if time bins are outside, we can't do anything, it's a property of the survey rather than the source
                mask = (idx_i==len(val['log_dt_bin_edges'])) + (idx_i==0)
                idx = idx_i[np.where((idx_i!=0)*(idx_i!=len(val['log_dt_bin_edges'])))]
                Y = Y[np.where((idx_i!=0)*(idx_i!=len(val['log_dt_bin_edges'])))]

            if 0 not in idx_i and (len(val['log_dt_bin_edges']) not in idx_i):
                idx = idx_i

            uniq_idx = np.unique(idx) 
            
            for k in range(len(uniq_idx)):
                mask = idx==uniq_idx[k]
                if val['log_dt_bin_edges'][uniq_idx[k]]<1.0: #dt less than 10 days
                    zero_dm_log = -2000.0 #to catch rapid rising objects (tuned value for ZTF)
                else:
                    zero_dm_log = -30.0
                    
                dm = Y[mask]
                indices = np.digitize(dm, val[str(uniq_idx[k])]['dm_bin_edges'])
                expectation = expected_D[key][str(uniq_idx[k])]
                if ((0 in indices) or (len(val[str(uniq_idx[k])]['dm_bin_edges']) in indices)):
                    mask = (indices==0) + (indices==len(val[str(uniq_idx[k])]['dm_bin_edges']))
                    MSK = indices[np.where((indices!=0)*(indices!=len(val[str(uniq_idx[k])]['dm_bin_edges'])))]
                    Sum = Sum + (np.sum(mask)) * (zero_dm_log-np.log10(expectation))
                    
                    if len(MSK)==0:
                        continue
                    

                if (0 not in indices) and (len(val[str(uniq_idx[k])]['dm_bin_edges']) not in indices):
                    MSK = indices

                try:
                    Sum = Sum + np.sum(np.log10(val[str(uniq_idx[k])]['dm_hist'][MSK-1]) - np.log10(expectation))
                except Exception as e:
                    print ("EXCEPTIONS", e)
                    print (f"The length of dm_hist is {len(val[str(uniq_idx[k])]['dm_hist'])}, and remaining indices to use is {MSK-1}")
                    
        score[K] = Sum #Score for each lightcurve id
            
    return score        

#=========================================================================================
      

        
        
        
def main(api_key, api_secret, antares_topic, passband=None, cc=None, thres=-3500.0):
    """
    Parameters
    ----------
    api_key : str
        Part of the API credentials for connecting to the Antares output Kafka stream.   
    api_secret : str
        Part of the API credentials for connecting to the Antares output Kafka stream.
    antares_topic : list of str
        The topics in the ANTARES output Kafka stream to be processed for anomalies.  
    passband : list
        The unique passbands of the light curve.  
    cc : dict
        The cross-passbands to be considered for computing the color-based dm-dt values.
        Example for ZTF is '{1:g-R, 2:R-g}'        
    thres : float
        The threshold to flag anomalous sources. Scores less than this value indicate anomalies.
        Note that this value is tuned using ZTF stellar variable sources. It may require tuning depending on your target of interest. 
    """            
    
    # ********** client declaration ********************    
    TOPICS = antares_topic
    CONFIG = {
    "api_key": api_key,
    "api_secret": api_secret,
    "group": 'ms_testanomalies', 
    "enable_auto_commit": False,
    }
    client = StreamingClient(TOPICS, **CONFIG)

    
    #### collect all (cross-)passbands to be considered
    grand_PB = []
    if passband is not None:
        grand_PB = grand_PB + passband
    if cc is not None:
        for K, V in cc.items():
            #grand_PB.append(V[0]+"-"+V[1])
            grand_PB.append(V)


    pdf_model = {}
    ##### read in the PDF model files    
    for K in grand_PB:
        pdf_model[K] = {}
        with open(f"{K}_logdt_dm.json", "rt") as f:
            temp = json.load(f)
        for key, Val in temp.items():
            if key=="log_dt_bin_edges":
                pdf_model[K][key] = np.array(Val)
            else:
                pdf_model[K][key] = {'dm_hist':np.array(Val['dm_hist']), 'dm_bin_edges':np.array(Val['dm_bin_edges'])}


    #### compute the normalization factor to be used in computing score
    expected_density = {}
    for key, val in pdf_model.items():
        expected_density[key] = {}
        for K, V in val.items():
            if K=='log_dt_bin_edges':
                continue
            dm_hist = np.array(V['dm_hist'])
            dm_bin_edges = np.array(V['dm_bin_edges'])
            expectation = np.sum(dm_hist**2.0*(dm_bin_edges[1:]-dm_bin_edges[0:-1]))
            expected_density[key][K] = expectation
            
            
    # ******************** main processing ********************
    for topic, locus in client.iter():
        ant_id = locus.locus_id  
        locus_ra = locus.ra
        locus_dec = locus.dec
        ztf_object_id = locus.properties['ztf_object_id']  # ZTF Object ID
        newest_alert_id = int(locus.properties['newest_alert_id'].split(":")[1])
        TS = locus.timeseries

        dets = TS[TS['ant_survey'] == 1]
        if len(dets) < 2:
            continue

        passband = dets['ant_passband'].data
        MAG = dets['ztf_magpsf'].data
        DMAG = dets['ztf_sigmapsf'].data
        MJD = dets['ant_mjd'].data


        mask = (MAG < 22.0) * (DMAG < 0.2) #remove bad alerts. Even bad NAN points are removed here automatically

        if np.sum(mask) == 0:
            continue

        thisobj = dmdt_model.dmt(ant_id, MJD[mask], MAG[mask], DMAG[mask], passband[mask], cc)

        b = thisobj.nocorr_bimod_structure_funct()

        out = {}
        out[thisobj.objname] = {}

        for i in range(len(grand_PB)):
            out[thisobj.objname][grand_PB[i]] = {"x":[-9999], "y":[-9999]}
        for key, val in b.items():
            if (key in grand_PB):
                out[thisobj.objname][key]["x"] = out[thisobj.objname][key]["x"] + val['log_dt']
                out[thisobj.objname][key]["y"] = out[thisobj.objname][key]["y"] + val['dmag']


        SCORE = get_score_normalized(out, pdf_model, expected_density)
        current_det_mag = MAG[np.argmax(MJD)]

        if list(SCORE.values())[0] < thres:
            print (f"HIT for antares locus {ant_id} with anomaly score {list(SCORE.values())[0]}")
            with open(f"anomalies_{topic.split('.')[-1]}.csv", "a") as f:
                f.write(f"{ant_id}, {topic}, {list(SCORE.values())[0]}\n")
        client.commit()

        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument('--api_key', type=str, required=True, help='ANTARES API credentials key')
    parser.add_argument('--api_secret', type=str, required=True, help='ANTARES API credentials secret')
    parser.add_argument('--antares_topic', type=str, required=True, help='ANTARES Kafka output topic(s) to subscribe')
    parser.add_argument('--passband', type=str, help='Passbands (e.g., R, g for ZTF)', default=None)
    parser.add_argument('--color_code', type=str, help="Example for ZTF is '{\"1\":\"g-R\", \"2\":\"R-g\"}'. The single and double quotes are required.", default=None)
    parser.add_argument('--score_threshold', type=float, help='Threshold for flagging anomalies', default=-3500.0)

    args = parser.parse_args()

    #print (args.api_key, args.api_secret, list(args.antares_topic.replace(' ','').split(',')), args.passband, args.color_code, args.score_threshold)    
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


    main(args.api_key, args.api_secret, list(args.antares_topic.replace(' ','').split(',')), passband=passband, cc=color_code, thres=args.score_threshold)
