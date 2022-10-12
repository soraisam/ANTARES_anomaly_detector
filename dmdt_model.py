import numpy as np
import itertools
import pandas as pd
import json


zero_log_dt = -7.0 # used as a filler below for simultaneous measurements across passbands

class dmt:
    """The `nocorr_bimod_structure_funct` method of this class computes the change in magnitude (dm) for a given change in time (dt) of a light curve."""
    
    def __init__(self, obj, date, mag, mag_err, passband, color_code, for_gen_pdf=False): 
        """
        Parameters
        ----------
        obj : str
            The locus_id / lightcurve_id
        date : numpy.ndarray
            The time stamps (MJD) of the light curve data points
        mag : numpy.ndarray
            The magnitudes of the light curve data points
        mag_err : numpy.ndarray
            The magnitude errors of the light curve data points
        passband : numpy.ndarray
            The passbands of the light curve data points  
        color_code : dict of {int : list}
            The passbands to be considered for computing the color-based dm-dt distributions.
            For ZTF, this is {1:["g", "R"], 2:["R", "g"]}
        for_gen_pdf : bool
            The boolean flag controlling whether the computed dm-dt values are to be output to a json file.
            Set this to True when generating the probability distributions of dm-dt values. 
        """        
    
        self.objname = obj 
        self.date = date.astype(float)
        self.mag = mag.astype(float)
        self.mag_err = mag_err.astype(float)
        if passband is not None:
            self.passband = passband.astype(str)
            self.filters = (np.unique(passband)).astype(str)            
        else:
            self.passband = passband
            self.filters = None
        self.color_code = color_code
        self.for_gen_pdf = for_gen_pdf

   
    def nocorr_bimod_structure_funct(self):
        out = {}
        F_series = {}
        if self.filters is not None:
            for F in list(self.filters):
                if (np.sum(self.passband==F)>=1):
                    date_per_F = self.date[self.passband==F]
                    mag_per_F = self.mag[self.passband==F]
                    magerr_per_F = self.mag_err[self.passband==F]

                    indices = np.argsort(date_per_F)
                    dates = date_per_F[indices]
                    mags = mag_per_F[indices]
                    magerr = magerr_per_F[indices]

                    ## need to remove duplicate light curve data points: they sometimes creep in for ZTF data as artifacts of Difference Imaging
                    mask = np.ones(len(dates), dtype=bool)
                    for k in range(len(dates)):
                        if list(dates).count(dates[k])>1:
                            mask[k] = False

                    dates = dates[mask]
                    mags = mags[mask]
                    magerr = magerr[mask]

                    if (len(dates)>=2):
                        log_dt = np.log10(np.diff(dates))
                        delta_mag = np.diff(mags)
                        dmag = np.array(delta_mag)
                        mask = ~np.isinf(log_dt)
                        log_dt = log_dt[mask]
                        dmag = dmag[mask]
                        indices = np.argsort(log_dt)
                        log_dt = log_dt[indices]
                        dmag = dmag[indices]
                        if (len(log_dt)>0):
                            out[F] = {"dmag":list(dmag), "log_dt":list(log_dt)}
                        else:
                            print (f"No dm-dt for filter {F}")

                    F_series[F] = {'t':dates, 'm':mags, 'e_m':magerr}

        if self.color_code is not None:
            for key, val in self.color_code.items():
                avail_pbs = list(F_series.keys())
                if val[0] not in avail_pbs or val[1] not in avail_pbs:
                    continue 

                idxes = pd.DataFrame(list(itertools.product(np.arange(len(F_series[val[0]]['t'])),
                                                          np.arange(len(F_series[val[1]]['t'])))), 
                                     columns=['l1', 'l2'])
                cross_dt = F_series[val[0]]['t'][idxes['l1'][:]]-F_series[val[1]]['t'][idxes['l2'][:]]
                cross_dm = F_series[val[0]]['m'][idxes['l1'][:]]-F_series[val[1]]['m'][idxes['l2'][:]]

                mask = cross_dt<0
                nonneg_dt = cross_dt[~mask]
                nonneg_dm = cross_dm[~mask]

                temp_log_dt=np.log10(nonneg_dt)
                idx=np.argsort(temp_log_dt)
                c_log_dt=temp_log_dt[idx]
                dc=nonneg_dm[idx]
                mask=np.isinf(c_log_dt)
                c_log_dt[mask]=zero_log_dt

                out[val[0]+"-"+val[1]] = {"log_dt":list(c_log_dt), "dmag":list(dc)}
            
        if self.for_gen_pdf == True:
            with open(f"./new_jsons/{self.objname}.json", "w") as f:
                json.dump(out, f, sort_keys=True, indent=4)
            return
            
        return out
