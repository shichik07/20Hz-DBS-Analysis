# -*- coding: utf-8 -*-
"""
Created on Wed May 10 17:56:20 2023

@author: doex9445
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import hddm
import patsy
import numpy as np
print(hddm.__version__)
os.chdir(r"C:\Users\doex9445\Dateien\Julius\20Hz\Data\Extracted")
os.getcwd()


data = hddm.load_csv('GoNoGo.csv')
# because we perform at least two model analyses - the first thing we want to see is how the Go trials and NoGo Cued trials differ by stimulation condition. For this analysis we ignore the NoGotrials and only at the Go trials

data_Go = data
# create a new column from existing column
data_Go['StopTrl'] = data_Go["GoNoGo"].apply(lambda x: 1 if "Stop" in x else 0)
data_Go['change_con'] = data_Go["GoNoGo"].apply(lambda x: 1 if "NoGo" in x else 0)
# filter Stop Trials
data_Go = data_Go[(data_Go.StopTrl != 1)]
# and rename a couple of columns
data_Go = data_Go.rename(columns={"RT":"rt", "Part_nr":"subj_idx", "Correct_Response": "response"})
data_Go.head(10)

#mat = patsy.dmatrix("0 + Stim + change_con + Stim:change_con", data_Go)
#mat = patsy.dmatrix("0 + C(Stim_verb, Treatment('OFF')):C(change_con, Treatment(0))", data_Go)
#mat_view = np.asarray(mat)

md_srt_stim = hddm.HDDMRegressor(data_Go, 
                                 "a ~ 0 + C(Stim_verb, Treatment('OFF')):C(change_con, Treatment(0))", # threshold depends on Stimulation condition and whether Go trials were cued
                                p_outlier = 0.05) # we say that 5% of our trials are outlier trials)
# find a good starting value
md_srt_stim.find_starting_values()
#get samples and discard a couple as burn ins
md_srt_stim.sample(10000, burn = 2000, )
md_srt_stim.save("GoNoGo_threshold_byc")

stats = md_srt_stim.gen_stats()
summary_Cued = stats[stats.index.isin(["a_C(Stim_verb, Treatment('OFF'))[130Hz]:C(change_con, Treatment(0))[1]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[20Hz]:C(change_con, Treatment(0))[1]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[OFF]:C(change_con, Treatment(0))[1]", 
                                  'v', 
                                  't'])]

summary_NoNCued = stats[stats.index.isin(["a_C(Stim_verb, Treatment('OFF'))[130Hz]:C(change_con, Treatment(0))[0]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[20Hz]:C(change_con, Treatment(0))[0]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[OFF]:C(change_con, Treatment(0))[0]", 
                                  'v', 
                                  't'])]
md_srt_stim.plot_posterior_predictive(figsize=(14,10))

Hz130, Hz20, OFF = md_srt_stim.nodes_db.loc[["a_C(Stim_verb, Treatment('OFF'))[130Hz]:C(change_con, Treatment(0))[1]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[20Hz]:C(change_con, Treatment(0))[1]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[OFF]:C(change_con, Treatment(0))[1]"], 'node']
hddm.analyze.plot_posterior_nodes([Hz130, Hz20, OFF])
plt.xlabel('decision threshold')
plt.ylabel('Posterior probability')
plt.title('Group mean posteriors of within-subject threshold effect Cued')
plt.legend(['130Hz', '20Hz', 'OFF'], 
           loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('hddm_threshold_GoNoGO_Cued.pdf')

print("P(20Hz > 130Hz Cued) = ", ((Hz20.trace() < Hz130.trace())).mean())
print("P(20Hz > OFF Cued) = ", ((Hz20.trace() <  OFF.trace())).mean())
print("P(OFF > 130Hz Cued) = ", ((OFF.trace() <  Hz130.trace())).mean())


Hz1301, Hz201, OFF1 = md_srt_stim.nodes_db.loc[["a_C(Stim_verb, Treatment('OFF'))[130Hz]:C(change_con, Treatment(0))[0]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[20Hz]:C(change_con, Treatment(0))[0]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[OFF]:C(change_con, Treatment(0))[0]"], 'node']
hddm.analyze.plot_posterior_nodes([Hz1301, Hz201, OFF1])
plt.xlabel('decision threshold')
plt.ylabel('Posterior probability')
plt.title('Group mean posteriors of within-subject threshold effect Not Cued')
plt.legend(['130Hz', '20Hz', 'OFF'], 
           loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('hddm_threshold_GoNoGO_NoNCued.pdf')

print("P(20Hz > 130Hz Non-Cued) = ", ((Hz201.trace() < Hz1301.trace())).mean())
print("P(20Hz > OFF Non-Cued) = ", ((Hz201.trace() <  OFF1.trace())).mean())
print("P(OFF > 130Hz Non-Cued) = ", ((OFF1.trace() <  Hz1301.trace())).mean())


ad = Hz20.trace()
ac = OFF.trace()
plt.hist(ad, density=True, bins = 90)
plt.hist(ac, density=True, bins = 90)

ab = Hz20.trace() - OFF.trace()
plt.hist(ab, density=True, bins = 90)
