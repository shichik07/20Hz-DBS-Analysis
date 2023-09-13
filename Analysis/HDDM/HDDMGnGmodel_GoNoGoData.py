# -*- coding: utf-8 -*-
"""
Created on Wed May 17 18:49:20 2023

@author: Julius Kricheldorff
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import hddm
from patsy import dmatrix
import numpy as np
import math
from joblib import Parallel, delayed
import random
import pickle
print(hddm.__version__)
os.chdir(r"C:\Users\doex9445\Dateien\Julius\20Hz\Data\Extracted")
os.getcwd()


data = hddm.load_csv('GoNoGo.csv')
# because we perform at least two model analyses - the first thing we want to see is how the Go trials and NoGo Cued trials differ by stimulation condition. For this analysis we ignore the NoGotrials and only at the Go trials
np.unique(data['Part_nr'])

# to see if the model works in principle we fit it to the particiapnts were we have equal amount of trials per condition first
#data2 = data[(data['Part_nr'] != 10)]

data_Ana = data[["RT", "Part_nr", "Stim_verb", "GoNoGo", "Correct_Response"]].copy()
# Then we need another column that shows us the correct response - 0 for nogo and 1 for Go trials - needed for model fitting stim_col parameter
data_Ana['stimulus'] = data_Ana["GoNoGo"].apply(lambda x: 0 if "Stop" in x else 1) # column to code Stop and Go trials
# Then we need another column that shows us the NoGo trial condition (including the Go trials)
data_Ana['ChangeCondition'] = data_Ana["GoNoGo"].apply(lambda x: 1 if "NoGo" in x else 0)

# trials of 0 second rts will be coded as NaN - or a NoGo Trial so I can set the next variable -awkward but works
data_Ana['RT_fake'] = data_Ana.RT.copy()
data_Ana['RT_fake'].loc[data_Ana['RT_fake'] == 0] = np.nan
# We need a response column that codes Stop responses as 0 and Go responses as 1
#data_Ana['response'] = data_Ana['RT_fake'].apply(lambda x: 0 if math.isnan(x)  else 1)
data_Ana['response'] = data_Ana.Correct_Response.copy() # data are accuracy coded, we using stimulus coding for the response variable

# and rename a couple of columns
data_Ana = data_Ana.rename(columns={"RT":"rt", "Part_nr":"subj_idx", "Correct_Response": "correct"})
#data_Go = data_Go[['rt', 'stim', 'response', 'Stim_verb']]
data_Ana.head(10)

# reset index, otherwise the link functions are a MESS just to be safe
data_Ana = data_Ana.reset_index()

# first the link function for the starting point bias - which we took from the hddm website
def z_link_func(x, data=data_Ana):
    stim = (np.asarray(dmatrix('0 + C(s, [[0], [1]])',
                              {'s': data.stimulus.loc[x.index]},return_type='dataframe'))
    )
    # Apply z = (1 - x) to flip them along 0.5
    z_flip = np.subtract(stim, x.to_frame())
    # The above inverts those values we do not want to flip,
    # so invert them back
    z_flip[stim == 0] *= -1
    return z_flip
    
# Second the Link function for the drift rate bias - which should be the easiest

def v_link_func(x, data=data_Ana):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.loc[x.index]})))
    return x * stim

# Okay so it might be possible that we do not need a separate link function for the decision threshold
# Instead, what we are going to do is create a new variable that contains Go trials (Go) versus NoGo (NoGo + Go) trials
# and estimate the threshold of the interaction with the normal stimulus function. I do not think the threshold for Go versus NoGo should vary
# Only by the condition. So lets set this up and see if we can recover our simulated parameter
v_reg = {'model': 'v ~ 1 + C(Stim_verb)', 'link_func': v_link_func}
z_reg = {'model': 'z ~ 1 + C(Stim_verb)', 'link_func': z_link_func}

a_reg = {'model': "a ~ 1 + C(Stim_verb, Treatment('OFF')):C(ChangeCondition, Treatment(0))", 'link_func': lambda x: x}



reg_model =  [z_reg, a_reg]
mod_zBias_threshold = hddm.HDDMRegressor(data_Ana, 
                             reg_model, 
                             include=['v', 'a', 'z', 't'])#,
                             #p_outlier = 0.05) # we say that 5% of our trials are outlier trials))

# find a good starting value
mod_zBias_threshold.find_starting_values()
#get samples and discard a couple as burn ins
mod_zBias_threshold.mcmc(dbname='GoNoGo_threshold_byc_Data2.db',db='pickle')
mod_zBias_threshold.sample(1000, burn = 50)

#get and save traces
traces = mod_zBias_threshold.nodes_db

model_samples = hddm.load('GoNoGo_threshold_byc_Data2.db')




staty = mod_zBias_threshold.print_stats()

Intercept_aOFF, NoGo_change, Hz130Go, Hz20Go, Hz130NoGo, Hz20NoGO = mod_zBias_threshold.nodes_db.loc[["a_Intercept",
                                            "a_C(ChangeCondition, Treatment(0))[T.1]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[T.130Hz]:C(ChangeCondition, Treatment(0))[0]", 
                                  "a_C(Stim_verb, Treatment('OFF'))[T.20Hz]:C(ChangeCondition, Treatment(0))[0]",
                                  "a_C(Stim_verb, Treatment('OFF'))[T.130Hz]:C(ChangeCondition, Treatment(0))[1]",
                                  "a_C(Stim_verb, Treatment('OFF'))[T.20Hz]:C(ChangeCondition, Treatment(0))[1]"
                                  ], 
                                           'node']
# Extract all other relevant parameter
Intercept_z130, z20, zOFF, v, t = mod_zBias_threshold.nodes_db.loc[["z_Intercept",
                                            "z_C(Stim_verb)[T.20Hz]", 
                                  "z_C(Stim_verb)[T.OFF]",
                                  "v",
                                  "t"], 
                                           'node']

# Next we construct the threshold off these parameters for each condition by stimulation
OffNoGo = Intercept_aOFF.trace() + NoGo_change.trace()
OffGo = Intercept_aOFF.trace() 
hz130NoGo = Intercept_aOFF.trace() + NoGo_change.trace() + Hz130NoGo.trace()
hz130Go = Intercept_aOFF.trace() + Hz130Go.trace()
Hz20NoGo = Intercept_aOFF.trace() + NoGo_change.trace() + Hz20NoGO.trace()
Hz20Go = Intercept_aOFF.trace() + Hz20Go.trace()

# same for the starting point bias
OFFZ = Intercept_z130.trace() + zOFF.trace()
Hz20z = Intercept_z130.trace() + z20.trace()

                                
# Plot and Store Distributions
distributions = [Hz20NoGo, Hz20Go, hz130NoGo,
                 hz130Go, OffNoGo,  OffGo]

# Plot densities using seaborn
sns.set(style="whitegrid")
for distribution in distributions:
    sns.kdeplot(distribution, shade=True)
    
plt.xlabel('decision threshold')
plt.ylabel('Posterior probability')
plt.title('Group mean posteriors of within-subject threshold effect Cued')
plt.legend(["Off change", "Off Go", '130Hz Change', "130Hz Go", '20Hz Change', '20Hz Go'], 
           loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('hddm_threshold_GoNoGO_Cued.pdf')

# Now lets perform posterior predictive checks by using the parameter estimates to simulate some data and see if we capture the characteristics of the data set
# code below is from the hddm wiki

def simulate_data_GnG(a, zz, v, t, z, sv=0, sz=0, st=0, condition = 0, nr_trials1=50, nr_trials2=25, nr_trials3=25):

    """
    Simulates stim-coded data.
    """

    parameters1 = {'a':a, 'v':v, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters2 = {'a':a + zz, 'v':v, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters3 = {'a':a + zz, 'v':v, 't':t, 'z':1-z, 'sv':sv, 'sz': sz, 'st': st}
    df_sim1, params_sim1 = hddm.generate.gen_rand_data(params=parameters1, size=nr_trials1, subjs=1, subj_noise=0)
    df_sim1['condition'] = "Go"
    df_sim2, params_sim2 = hddm.generate.gen_rand_data(params=parameters2, size=nr_trials2, subjs=1, subj_noise=0)
    df_sim2['condition'] = "Go - NoGo"
    df_sim3, params_sim3 = hddm.generate.gen_rand_data(params=parameters3, size=nr_trials3, subjs=1, subj_noise=0)
    df_sim3['condition'] = "NoGo"
    df_sim = pd.concat((df_sim1, df_sim2, df_sim3))
    df_sim['bias_response'] = df_sim.apply(get_choice, 1)
    df_sim['correct'] = df_sim['response'].astype(int)
    df_sim['response'] = df_sim['bias_response'].astype(int)
    df_sim['stimulus'] = np.array((np.array(df_sim['response']==1) & np.array(df_sim['correct']==1)) + (np.array(df_sim['response']==0) & np.array(df_sim['correct']==0)), dtype=int)
    df_sim['StimC'] = condition
    df_sim = df_sim.drop(columns=['bias_response'])

    return df_sim

def get_choice(row):

    #if row.condition == 'present':
    if row.condition == 'Go' or row.condition == "Go - NoGo":
        if row.response == 1:
            return 1
        else:
            return 0
    #elif row.condition == 'absent':
    elif row.condition == "NoGo":
        if row.response == 0:
            return 1
        else:
            return 0
        
# simulate data
random.seed(852)
# settings
go_nogo = True # should we put all RTs for one choice alternative to NaN (go-no data)?
n_subjects = 16
GoTrials = 75
NoGoTrials =25

# parameters:
params0 = {'cond':'OFF', 'v': v.trace().mean(), 'a': OffGo.mean(), 'zz': OffNoGo.mean() - OffGo.mean(),'t': t.trace().mean(), 'z': OFFZ.mean(), 'sz':0, 'st':0, 'sv':0}
params1 = {'cond':'130Hz', 'v': v.trace().mean(), 'a': hz130Go.mean(), 'zz': hz130NoGo.mean() - hz130Go.mean(),'t': t.trace().mean(), 'z': Intercept_z130.trace().mean(), 'sz':0, 'st':0, 'sv':0}
params2 = {'cond':'20Hz', 'v': v.trace().mean(), 'a': Hz20Go.mean(), 'zz': Hz20NoGo.mean() - Hz20Go.mean(),'t': t.trace().mean(), 'z': Hz20z.mean(), 'sz':0, 'st':0, 'sv':0}


# # parameters:
# params0 = {'cond':'OFF', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':0.7, 'dc':0, 'sz':0, 'st':0, 'sv':0}
# params1 = {'cond':'130Hz', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':0.7, 'dc':0, 'sz':0, 'st':0, 'sv':0}
# params2 = {'cond':'20Hz', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':.7, 'dc':0, 'sz':0, 'st':0, 'sv':0}


# simulate:
dfs = []
for i in range(n_subjects):
    df0 = simulate_data_GnG(z=params0['z'], a=params0['a'], v=params0['v'],
                        t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'],
                        condition=params0['cond'], zz = params0['zz'])
    df1 = simulate_data_GnG(z=params1['z'], a=params1['a'], v=params1['v'], 
                        t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'],
                        condition=params1['cond'], zz = params1['zz'])
    df2 = simulate_data_GnG(z=params2['z'], a=params2['a'], v=params2['v'],
                        t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'],
                        condition=params2['cond'], zz = params2['zz'])
    df = pd.concat((df0, df1, df2))
    df['subj_idx'] = i
    dfs.append(df)

# combine in one dataframe:
df_emp = pd.concat(dfs)
go_nogo = 1
if go_nogo:
    df_emp.loc[df_emp["response"]==0, 'rt'] = np.NAN

# GoNoGo_zBias_threshold = hddm.HDDMRegressor(data_Go, 
#                                  "a ~ 0 + C(Stim_verb, Treatment('OFF')):C(change_con, Treatment(0))", # threshold depends on Stimulation condition and whether Go trials were cued
#                                 p_outlier = 0.05) # we say that 5% of our trials are outlier trials)
# # find a good starting value
# md_srt_stim.find_starting_values()
# #get samples and discard a couple as burn ins
# md_srt_stim.sample(10000, burn = 2000, )
# md_srt_stim.save("GoNoGo_threshold_byc")

# def fit_subject(data, quantiles):

#     """
#     Simulates stim-coded data.
#     """

#     subj_idx = np.unique(data['subj_idx'])
#     m = hddm.HDDMStimCoding(data_Go, stim_col='stim', split_param='v', drift_criterion=True, bias=True, p_outlier=0,
#                             depends_on={'v':'Stim_verb', 'a':'Stim_verb', 't':'Stim_verb', 'z':'Stim_verb', 'dc':'Stim_verb', })
#     m.optimize('gsquare', quantiles=quantiles, n_runs=8)
#     res = pd.concat((pd.DataFrame([m.values], index=[subj_idx]), pd.DataFrame([m.bic_info], index=[subj_idx])), axis=1)
#     return res

# n_subjects = len(np.unique(data_Go['subj_idx']))

# params_fitted = pd.concat(Parallel(n_jobs=4)(delayed(fit_subject)(data[1], quantiles)
#                                                       for data in data_Go.groupby('subj_idx')))


