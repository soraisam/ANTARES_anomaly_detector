**Overview**

This codebase (compatible with python 3.x) contains the production version of the anomaly detection algorithm presented by <a href="https://ui.adsabs.harvard.edu/abs/2020ApJ...892..112S/abstract">Soraisam et al. 2020</a> for processing ZTF alerts. It is running on the alert-broker <a href="https://antares.noirlab.edu">Antares</a>. 

If you make use of the codes here, please cite the above paper (including for modified versions of the codes).  

**Contents**

- `anom_prod.py`: main source code for computing the anomaly score of an incoming alert.

- `dmdt_model.py`: source code for computing the dm-dt of a given light curve.

- `*.json`: probability distributions used for computing the anomaly scores of incoming alerts. An example file `g_logdt_dm.json` is included. Make sure to generate similar distributions appropriate for your application.

- `gen_pdfs.py`: source code to generate the above probability distributions.

- `new_jsons/`: folder where the dm-dt values for a given light curve in the *reference* sample will be written.

**Workflow**

1. Compile a sizable *reference* sample of loci (i.e., light curves) against which the peculiarity of a given real-time alert is to be assessed. 

2. Use the `dmdt_model` module to compute the dm-dt json file for each light curve of the reference sample. The output file will be written in `new_jsons/`.<br>
    <ins>Example</ins><br>
    ```
    from dmdt_model import dmt
    import pandas as pd 
    
    LC = pd.read_csv("example_lc_48904.csv")
    p = dmt(obj='example_lc_48904', date=LC['mjd'].values, mag=LC['mag'].values, mag_err=LC['magerr'].values,
      	 passband=LC['pb'].values, color_code={1:["g", "R"], 2:["R", "g"]}, for_gen_pdf=True)
    p.nocorr_bimod_structure_funct()
    ```
3. Generate conditional probability distributions P(dm|dt) for the reference sample using the `gen_pdfs.py` script. This is typically a one-off process provided the reference sample does not considerably evolve. If one has a reason to believe otherwise, regenerate the distributions at an appropriate frequency (e.g., every few months, every few weeks, etc.).<br>
    <ins>Example</ins><br>
    ```
    python gen_pdfs.py --passband 'g R' --color_code '{"1":"g-R", "2":"R-g"}'
    ```
    The above will generate the distributions for individual passbands g and R as well as for cross-passbands (or pseudo-color) g-R and R-g. 
    
4. Finally, run the script `anom_prod.py` (for example as a cron job). Note that this is currently configured for running downstream of Antares, processing one (or more) of its output alert streams (i.e., Antares Kafka topics).<br>
    <ins>Example</ins><br>
    ```
    python anom_prod.py --help # for examining the usage (see below)
    python anom_prod.py --api_key test --api_secret blahXXXXTTTTTT --antares_topic high_snr_staging --passband 'gR' --color_code '{"1":"g-R", "2":"R-g"}' --score_threshold -4000.0
    ```
    The threshold used to flag anomalous sources can be tuned. Use the output file `anomalies_<topic>.csv` whose last column contains the scores to tune the `score_threshold` if needed.   

