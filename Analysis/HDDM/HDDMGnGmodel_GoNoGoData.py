# -*- coding: utf-8 -*-
"""
Created on Wed May 17 18:49:20 2023

@author: Julius Kricheldorff
"""

# git repo adress C/Users/doex9445/Dateien/Julius/20Hz-DBS-Analysis/Analysis

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
os.chdir(r"C:\Users\doex9445\Dateien\Julius\20Hz-DBS-Analysis\Data\Extracted")
os.getcwd()

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

def simulate_data_GnG(a, vNoGoGo, vStop, vGo, t, z, sv=0, sz=0, st=0, condition = 0, nr_trials1=100, nr_trials2=50, nr_trials3=50):

    """
    Simulates stim-coded data.
    """

    parameters1 = {'a':a, 'v': vGo, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters2 = {'a':a, 'v': vNoGoGo, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st} # we set the parameter to 0.5 here because we expect go and nogo trials here to be equally likely
    parameters3 = {'a':a, 'v': vStop, 't':t, 'z':1-z, 'sv':sv, 'sz': sz, 'st': st}
    df_sim1, params_sim1 = hddm.generate.gen_rand_data(params=parameters1, size=nr_trials1, subjs=1, subj_noise=0)
    df_sim1['condition'] = "Go"
    df_sim2, params_sim2 = hddm.generate.gen_rand_data(params=parameters2, size=nr_trials2, subjs=1, subj_noise=0)
    df_sim2['condition'] = "Go - NoGo"
    df_sim3, params_sim3 = hddm.generate.gen_rand_data(params=parameters3, size=nr_trials3, subjs=1, subj_noise=0)
    df_sim3['condition'] = "NoGo"
    df_sim = pd.concat((df_sim1, df_sim2, df_sim3))
    df_sim['correct'] = df_sim.apply(get_choice, 1)
    df_sim['StimC'] = condition
    
    return df_sim

# function to plot the data
def seaborn_data_plotting(data):
    sns.set_style("whitegrid")
    pal = sns.color_palette('Set2')
    #set figure size
    sns.set(rc = {'figure.figsize': (10, 6)})
    #data = pd.DataFrame(df_emp)
    data.condition.astype('category')
    #reset index
    data = data.reset_index()
    data["error"] = 1 - data["correct"]
    #set order so the medians match
    hue_ord = ['130Hz', '20Hz', 'OFF']
    # set up subplots
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    # plot RT
    box_plot = sns.boxplot(data = data[data['condition'] != 'NoGo - Stop'],
                           hue ='condition', y = 'rt', x = 'StimC', 
                           order=hue_ord, ax = ax1, palette= pal)
    #sns.despine(offset = 10, trim = True) # make the axis more exiting
    # median RT added to plot
    median_RT = data[data['condition'] != 'NoGo - Stop'].groupby(['StimC', 'condition'])['rt'].median().round(3)
    vertical_offset = data[data['condition'] != 'NoGo - Stop'].rt.median()* 0.05 # offset from median for display
    tick_range = [0,0,1,1,2,2] # ticks
    counter = 0
    
    # write median values in graph
    for xtick in tick_range:
        if counter%2 ==1:
            offset = 0.2
        else: 
            offset = -0.2
        
        box_plot.text(xtick + offset, median_RT[counter] + vertical_offset, median_RT[counter], 
                      horizontalalignment='center',size='x-small',color='w',weight='bold', bbox=dict(facecolor='#445A64'))
        counter +=1
        
    # plot accuracy values
    ax2 = fig.add_subplot(122)
    mean_err = data.groupby(['StimC', 'condition'])['error'].mean().round(3)
    bar_plot = sns.barplot(data = data, hue = 'condition', y = 'error',  
                x = 'StimC',  order=hue_ord,ax = ax2, palette= pal)
    
    vertical_offset = data.error.mean()* 0.05 # offset from median for display
    tick_range_err = [0,0,0,1,1,1,2,2,2] # ticks
    counter = 0

    # write median values in graph
    for xtick in tick_range_err:
        if counter%3 ==1:
            offset = 0
        elif counter%3 == 2:
            offset = 0.25
        else: 
            offset = -0.25
        
        bar_plot.text(xtick + offset, mean_err[counter] + vertical_offset, mean_err[counter], 
                      horizontalalignment='center',size='x-small',color='w',weight='bold', bbox=dict(facecolor='#445A64'))
        counter +=1
    
    plt.close(2)
    plt.close(3)
    plt.tight_layout()

data = hddm.load_csv('GoNoGo.csv')
# because we perform at least two model analyses - the first thing we want to see is how the Go trials and NoGo Cued trials differ by stimulation condition. For this analysis we ignore the NoGotrials and only at the Go trials
np.unique(data['Part_nr'])

# to see if the model works in principle we fit it to the particiapnts were we have equal amount of trials per condition first
#data2 = data[(data['Part_nr'] != 10)]

# only use columns we actually need
data_Ana = data[["RT", "Part_nr", "Stim_verb", "GoNoGo", "Correct_Response"]].copy()
# filter reaction times larger than 3s as we did in the main analysis
data_Ana = data_Ana[data_Ana["RT"] < 3]
# Then we need another column that shows us the correct response - 0 for nogo and 1 for Go trials - needed for model fitting stim_col parameter
data_Ana['stimulus'] = data_Ana["GoNoGo"].apply(lambda x: "NoGo" if "Stop" in x else "Go") # column to code Stop and Go trials

# just to make things clear on the RT side of things
data_Ana.RT[data.RT == 0] = -1
# get a response column for wheather participants responded Go or NoGo
data_Ana['response'] = data_Ana["RT"].apply(lambda x: 0 if x==-1.0 else 1) # column to code Stop and Go trials
# Then we need another column that shows us the NoGo trial condition (including the Go trials)
data_Ana['ChangeCondition'] = data_Ana["GoNoGo"].apply(lambda x: 1 if "NoGo" in x else 0)
# get a sum contrast for the Go -NoGo vs Go condition 
data_Ana['Contrast_Go'] = data_Ana["GoNoGo"].apply(lambda x: 1 if x == "Go" else(-1 if x == "NoGo - Go" else 0))
# Then we need another column that shows us the NoGo trial condition (including the Go trials)
# get a sum contrast for the NoGo vs Go condition 
data_Ana['Contrast_GovsStop'] = data_Ana["GoNoGo"].apply(lambda x: -1 if x == "NoGo - Stop" else 1)
 

# and rename a couple of columns
data_Ana = data_Ana.rename(columns={"RT":"rt", "Part_nr":"subj_idx", "Correct_Response": "correct", "Stim_verb": "StimC", "GoNoGo":"condition"})
# reset index, otherwise the link functions are a MESS just to be safe
data_Ana = data_Ana.reset_index(drop = True)
#look at data
data_Ana.head(10)
# plot data
seaborn_data_plotting(data_Ana)
# reset index, otherwise the link functions are a MESS just to be safe
data_Ana = data_Ana.reset_index(drop = True)

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

v_reg = {'model': 'v ~ 1 + Contrast_GovsStop', 'link_func': lambda v: v}
z_reg = {'model': 'z ~ 1 + Contrast_Go', 'link_func': z_link_func}
a_reg = {'model': "a ~ 0 + C(StimC)", 'link_func': lambda a: a}

reg_model =  [z_reg, v_reg, a_reg]
model_driftbias = hddm.HDDMRegressor(data_Ana, 
                             reg_model, 
                             include=['z' , 'v'])

#get samples and discard a couple as burn ins
model_driftbias.sample(500, burn = 100, dbname='GoNoGo_threshold_byc_Sim',db='pickle')
# save model
model = hddm.load('GoNoGo_threshold_byc_Sim')
#check stats
staty = model_driftbias.print_stats()
#get samples and discard a couple as burn ins
#mod_zBias_threshold.mcmc(dbname='GoNoGo_threshold_byc_Data2.db',db='pickle')
#mod_zBias_threshold.sample(1000, burn = 50)


t, z, v_Go, v_NoGoGo, v_Stop , a_130, a_20, a_OFF = model_driftbias.nodes_db.loc[["t", "z_Intercept", "v_Intercept", "v_C(condition)[T.NoGo - Go]",
                                                         "v_C(condition)[T.NoGo - Stop]", "a_Intercept", "a_C(StimC)[T.20Hz]", "a_C(StimC)[T.OFF]"], 
                                                                  'node']
# use the mean as the parameter of interest
par_t = t.trace()
par_z = 1 - z.trace() # because stop was used as the reference
par_vGo = v_Go.trace()
par_vNoGoGo = v_NoGoGo.trace() + par_vGo
par_vStop = v_Stop.trace() + par_vGo
par_a130 = a_130.trace()
par_a20 = a_20.trace() + par_a130
par_aOFF =  a_OFF.trace() + par_a130

# next we check if the model parameter fit the data

# parameters:
n_subjects = 16
params0 = {'cond':'OFF', 'vGo': par_vGo.mean(), 'vNoGoGo': par_vNoGoGo.mean(), 'vStop': par_vStop.mean(), 
           'a':par_aOFF.mean(), 't':par_t.mean(), 'z':par_z.mean(), 'sz':0, 'st':0, 'sv':0}
params1 = {'cond':'130Hz', 'vGo': par_vGo.mean(), 'vNoGoGo': par_vNoGoGo.mean(), 'vStop': par_vStop.mean(), 
           'a':par_a130.mean(), 't':par_t.mean(), 'z':par_z.mean(), 'sz':0, 'st':0, 'sv':0}
params2 = {'cond':'20Hz', 'vGo': par_vGo.mean(), 'vNoGoGo': par_vNoGoGo.mean(), 'vStop': par_vStop.mean(), 
           'a':par_a20.mean(), 't':par_t.mean(), 'z':par_z.mean(), 'sz':0, 'st':0, 'sv':0}

# simulate complex task:
dfs = []
for i in range(n_subjects):
    df0 = simulate_data_GnG(z=params0['z'], a=params0['a'], vGo=params0['vGo'], vNoGoGo = params0['vNoGoGo'], vStop = params0['vStop'],
                        t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'], 
                        condition=params0['cond'])
    df1 = simulate_data_GnG(z=params1['z'], a=params1['a'], vGo=params1['vGo'], vNoGoGo = params1['vNoGoGo'], vStop = params1['vStop'],
                        t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'],
                        condition=params1['cond'])
    df2 = simulate_data_GnG(z=params2['z'], a=params2['a'], vGo=params2['vGo'], vNoGoGo = params2['vNoGoGo'], vStop = params2['vStop'],
                        t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'],
                        condition=params2['cond'])
    df = pd.concat((df0, df1, df2))
    df['subj_idx'] = i
    dfs.append(df)

# combine in one dataframe:
df_emp = pd.concat(dfs)
df_emp.loc[df_emp["response"]==0, 'rt'] = -1
    
# Plot the accuracy data for Go and NoGo responses
seaborn_data_plotting(df_emp)


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


