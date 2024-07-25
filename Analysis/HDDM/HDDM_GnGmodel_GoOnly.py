# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 14:53:49 2023

@author: doex9445
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
        
def posterior_sim(model, seed, nr_participants = 17, runs = 500):
    """
    Draw from the posterior and simulate data
    
    Args: 
        model: hddm model
        nr_participants: integer - dicating for how many participants we simulate
        runs: integer - how many draws from the posterior we use 
        
    """

    # pull parameter estimates
    t, v, a_130, a_Go_NoGo, a_OFF_Go, a_OFF_NoGo, a_20_Go, a_20_NoGo, v_std, t_std, a_std = model.nodes_db.loc[["t", "v", "a_Intercept",
                                                              "a_C(condition)[T.NoGo - Go]","a_C(StimC)[T.OFF]:C(condition)[Go]", 
                                                              "a_C(StimC)[T.OFF]:C(condition)[NoGo - Go]","a_C(StimC)[T.20Hz]:C(condition)[Go]",
                                                              "a_C(StimC)[T.20Hz]:C(condition)[NoGo - Go]", "v_std", "t_std", "a_Intercept_std" 
                                                              ], 
                                                                      'node']
    # get number of posterior draws
    n_draws = len(t.trace())
    
    # set seed so this can be reproduced
    random.seed(seed)
    
    # get 500 samples
    samples = random.sample(range(n_draws), runs)
    
    
    # get trace of the parameter
    par_t = t.trace()
    par_v = v.trace()
    par_a130 = a_130.trace()
    par_aNoGoGo = a_Go_NoGo.trace() 
    par_aOFF_Go = par_a130 + a_OFF_Go.trace()
    par_aOFF_NoGo = par_a130 + par_aNoGoGo + a_OFF_NoGo.trace()
    par_a20_Go = par_a130 + a_20_Go.trace()
    par_a20_NoGo = par_a130 + par_aNoGoGo + a_20_NoGo.trace()
    par_a130_Go = par_a130
    par_a130_NoGo = par_a130 + par_aNoGoGo
    par_sv = v_std.trace()
    par_st = t_std.trace()
    par_sa = a_std.trace()
    
    # create an empty dataframe where we save our variables
    columns = ["rt", "error", "condition", "draw", "StimC"]
    emp_df = pd.DataFrame(columns = columns)
    
    # run a for loop for your nr of runs
    for i in range(runs):
        smp = samples[i]

        # parameters:
        n_subjects = nr_participants
        params0 = {'cond':'OFF', 'aGo': par_aOFF_Go[smp], 'aNoGoGo': par_aOFF_NoGo[smp],
                   'v':par_v[smp], 't':par_t[smp], 'sz':0, 'st':par_st[smp], 'sv': par_sv[smp], 'sa': par_sa[smp]}
        params1 = {'cond':'130Hz', 'aGo': par_a130_Go[smp], 'aNoGoGo': par_a130_NoGo[smp],
                   'v':par_v[smp], 't':par_t[smp], 'sz':0, 'st':par_st[smp], 'sv': par_sv[smp], 'sa': par_sa[smp]}
        params2 = {'cond':'20Hz', 'aGo': par_a20_Go[smp], 'aNoGoGo': par_a20_NoGo[smp],
                   'v':par_v[smp], 't':par_t[smp], 'sz':0, 'st':par_st[smp], 'sv': par_sv[smp], 'sa': par_sa[smp]}
    
    
        # simulate complex task:
        dfs = []
        for i in range(n_subjects):
            df0 = simulate_data_GnG(v=params0['v'], aGo=params0['aGo'],aNoGo=params0['aNoGoGo'],
                                t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'], sa=params0['sa'],
                                condition=params0['cond'])
            df1 = simulate_data_GnG(v=params1['v'], aGo=params1['aGo'],aNoGo=params1['aNoGoGo'],
                                t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'], sa=params0['sa'],
                                condition=params1['cond'])
            df2 = simulate_data_GnG(v=params2['v'], aGo=params2['aGo'],aNoGo=params2['aNoGoGo'],
                                t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'], sa=params0['sa'],
                                condition=params2['cond'])
            df = pd.concat((df0, df1, df2))
            df['subj_idx'] = i
            dfs.append(df)
    
        # combine in one dataframe:
        df_emp = pd.concat(dfs)
        df_emp["error"] = 1 - df_emp.correct
        
        # save values in dataframe
        rt_data = df_emp.groupby(["StimC", "condition"]).rt.mean().reset_index()
        # save values in dataframe
        error_data = df_emp.groupby(["StimC", "condition"]).error.mean().reset_index()
        # merge
        new_df = pd.merge(rt_data, error_data, on = ["StimC", "condition"])
        # add column with draw
        new_df["draw"] = smp
        
        # save final result
        emp_df = emp_df.append(new_df, ignore_index=True)
        
    return(emp_df)

def simulate_data_GnG(aGo, aNoGo, v, t, sa, sv=0, sz=0, st=0,  condition = 0, nr_trials1=75, nr_trials2=25, nr_trials3=25):

    """
    Simulates stim-coded data.
    """

    parameters1 = {'a':aGo, 'v': v, 't':t, 'z': 0.5, 'sv':sv, 'sz': sz, 'st': st, 'sa': sa}
    parameters2 = {'a':aNoGo, 'v': v, 't':t, 'z': 0.5, 'sv':sv, 'sz': sz, 'st': st, 'sa': sa} # we set the parameter to 0.5 here because we expect go and nogo trials here to be equally likely
    df_sim1, params_sim1 = hddm.generate.gen_rand_data(params=parameters1, size=nr_trials1, subjs=1, subj_noise=0)
    df_sim1['condition'] = "Go"
    df_sim2, params_sim2 = hddm.generate.gen_rand_data(params=parameters2, size=nr_trials2, subjs=1, subj_noise=0)
    df_sim2['condition'] = "Go - NoGo"
    df_sim = pd.concat((df_sim1, df_sim2))
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
    tick_range_err = [0,0,0,1,1,1] # ticks
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
data_Ana = data_Ana[(data_Ana["RT"] < 3) & (data_Ana["RT"] >= 0.2)]
# Then we need another column that shows us the correct response - 0 for nogo and 1 for Go trials - needed for model fitting stim_col parameter
#data_Ana['stimulus'] = data_Ana["GoNoGo"].apply(lambda x: "NoGo" if "Stop" in x else "Go") # column to code Stop and Go trials

# just to make things clear on the RT side of things
#data_Ana.RT[data.RT == 0] = -1
# get a response column for wheather participants responded Go or NoGo
data_Ana['response'] = data_Ana["Correct_Response"]
# Then we need another column that shows us the NoGo trial condition (including the Go trials)
data_Ana['ChangeCondition'] = data_Ana["GoNoGo"].apply(lambda x: 1 if "NoGo" in x else 0)
# get a sum contrast for the Go -NoGo vs Go condition 
data_Ana['Contrast_Go'] = data_Ana["GoNoGo"].apply(lambda x: 1 if x == "Go" else(-1 if x == "NoGo - Go" else 0))
# Then we need another column that shows us the NoGo trial condition (including the Go trials)
# get a sum contrast for the NoGo vs Go condition 
data_Ana['Contrast_GovsStop'] = data_Ana["GoNoGo"].apply(lambda x: -1 if x == "NoGo - Stop" else 1)
 

# and rename a couple of columns
data_Ana = data_Ana.rename(columns={"RT":"rt", "Part_nr":"subj_idx", "Stim_verb": "StimC","Correct_Response" : "correct", "GoNoGo":"condition"})
# reset index, otherwise the link functions are a MESS just to be safe
data_Ana = data_Ana.reset_index(drop = True)
#look at data
data_Ana.head(10)
# plot data
seaborn_data_plotting(data_Ana)
# Filter StopTrials out
data_Ana = data_Ana[data_Ana["condition"] != "NoGo - Stop"]

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

#v_reg = {'model': 'v ~ 1 + C(StimC):C(condition)', 'link_func': lambda v: v}
#z_reg = {'model': 'z ~ 1 + C(condition)', 'link_func':  lambda v: v}
a_reg = {'model': "a ~ 0 +  C(StimC):C(condition)", 'link_func': lambda a: a}

reg_model =  [a_reg]
model_driftbias = hddm.HDDMRegressor(data_Ana, 
                              reg_model, 
                              include=['v', 't'])

# model_driftbias = hddm.HDDM(data_Ana, bias=True, include=('v', 'a', 'z', 't', ))

# model_driftbias = hddm.HDDM(data_Ana, include=('v', 'a', 'z', 't', ), depends_on={'a': {'StimC', 'condition'}})


# model_driftbias = hddm.HDDMStimCoding(data_Ana,stim_col='condition', include=('z', 't'), 
#                         depends_on={'a': 'StimC', 'v': 'condition', 'z': 'condition'})

#get samples and discard a couple as burn ins
model_driftbias.sample(10000, burn = 2000, dbname='GoNoGo_threshold_byc_Sim',db='pickle')

model_driftbias.poster

#check stats
staty = model_driftbias.print_stats()

# pull parameter estimates
t, v, a_130, a_Go_NoGo, a_OFF_Go, a_OFF_NoGo, a_20_Go, a_20_NoGo, v_std, t_std, a_std = model_driftbias.nodes_db.loc[["t", "v", "a_Intercept",
                                                          "a_C(condition)[T.NoGo - Go]", "a_C(StimC)[T.OFF]:C(condition)[Go]", 
                                                          "a_C(StimC)[T.OFF]:C(condition)[NoGo - Go]","a_C(StimC)[T.20Hz]:C(condition)[Go]",
                                                          "a_C(StimC)[T.20Hz]:C(condition)[NoGo - Go]", "v_std", "t_std", "a_Intercept_std" 
                                                          ], 
                                                                  'node']
# use the mean as the parameter of interest
par_t = t.trace()
par_v = v.trace()
par_a130 = a_130.trace()
par_aNoGoGo = a_Go_NoGo.trace() 
par_aOFF_Go = par_a130 + a_OFF_Go.trace()
par_aOFF_NoGo = par_a130 + par_aNoGoGo + a_OFF_NoGo.trace()
par_a20_Go = par_a130 + a_20_Go.trace()
par_a20_NoGo = par_a130 + par_aNoGoGo + a_20_NoGo.trace()
par_a130_Go = par_a130
par_a130_NoGo = par_a130 + par_aNoGoGo
par_sv = v_std.trace()
par_st = t_std.trace()
par_sa = a_std.trace()

par_a20_NoGo.mean()
par_aOFF_NoGo.mean()

ab = par_a20_NoGo - par_aOFF_NoGo
plt.hist(ab, density=True, bins = 90)
# parameters:
n_subjects = 160
params0 = {'cond':'OFF', 'aGo': par_aOFF_Go.mean(), 'aNoGoGo': par_aOFF_NoGo.mean(),
           'v':par_v.mean(), 't':par_t.mean(), 'sz':0, 'st':par_st.mean(), 'sv': par_sv.mean(), 'sa': par_sa.mean()}
params1 = {'cond':'130Hz', 'aGo': par_a130_Go.mean(), 'aNoGoGo': par_a130_NoGo.mean(),
           'v':par_v.mean(), 't':par_t.mean(), 'sz':0, 'st':par_st.mean(), 'sv': par_sv.mean(), 'sa': par_sa.mean()}
params2 = {'cond':'20Hz', 'aGo': par_a20_Go.mean(), 'aNoGoGo': par_a20_NoGo.mean(),
           'v':par_v.mean(), 't':par_t.mean(), 'sz':0, 'st':par_st.mean(), 'sv': par_sv.mean(), 'sa': par_sa.mean()}


# simulate complex task:
dfs = []
for i in range(n_subjects):
    df0 = simulate_data_GnG(v=params0['v'], aGo=params0['aGo'],aNoGo=params0['aNoGoGo'],
                        t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'], sa=params0['sa'],
                        condition=params0['cond'])
    df1 = simulate_data_GnG(v=params1['v'], aGo=params1['aGo'],aNoGo=params1['aNoGoGo'],
                        t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'], sa=params0['sa'],
                        condition=params1['cond'])
    df2 = simulate_data_GnG(v=params2['v'], aGo=params2['aGo'],aNoGo=params2['aNoGoGo'],
                        t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'], sa=params0['sa'],
                        condition=params2['cond'])
    df = pd.concat((df0, df1, df2))
    df['subj_idx'] = i
    dfs.append(df)

# combine in one dataframe:
df_emp = pd.concat(dfs)
#df_emp.loc[df_emp["response"]==0, 'rt'] = -1

seaborn_data_plotting(df_emp)

posterior_pred = posterior_sim(model = model_driftbias , seed = 123, nr_participants = 17, runs = 500)

rt_data = posterior_pred.groupby(["StimC", "condition"]).rt.mean().reset_index()
# save values in dataframe
error_data = posterior_pred.groupby(["StimC", "condition"]).error.mean().reset_index()

# # pull parameter estimates
# t, z, z_NoGoGo, v_Go, v_NoGoGo,  a_130, a_20, a_OFF = model_driftbias.nodes_db.loc[["t", "z_Intercept", "z_C(condition)[T.NoGo - Go]", 
#                                                                                             "v_Intercept", "v_C(condition)[T.NoGo - Go]",
#                                                           "a_Intercept", "a_C(StimC)[T.20Hz]", "a_C(StimC)[T.OFF]"], 
#                                                                   'node']
# # use the mean as the parameter of interest
# par_t = t.trace()
# par_zGo = z.trace() 
# par_zNoGoGo = z.trace() + z_NoGoGo.trace()
# par_vGo = v_Go.trace()
# par_vNoGoGo = v_NoGoGo.trace() + par_vGo
# par_a130 = a_130.trace()
# par_a20 = a_20.trace() + par_a130
# par_aOFF =  a_OFF.trace() + par_a130
# simulate data

# parameters:
n_subjects = 16
params0 = {'cond':'OFF', 'vGo': par_vGo.mean(), 'vNoGoGo': par_vNoGoGo.mean(),
           'a':par_aOFF.mean(), 't':par_t.mean(), 'par_zGo':par_zGo.mean(), 'par_zNoGoGo':par_zNoGoGo.mean(), 'sz':0, 'st':0, 'sv':0}
params1 = {'cond':'130Hz', 'vGo': par_vGo.mean(), 'vNoGoGo': par_vNoGoGo.mean(), 
           'a':par_a130.mean(), 't':par_t.mean(), 'par_zGo':par_zGo.mean(), 'par_zNoGoGo':par_zNoGoGo.mean(),'sz':0, 'st':0, 'sv':0}
params2 = {'cond':'20Hz', 'vGo': par_vGo.mean(), 'vNoGoGo': par_vNoGoGo.mean(), 
           'a':par_a20.mean(), 't':par_t.mean(), 'par_zGo':par_zGo.mean(), 'par_zNoGoGo':par_zNoGoGo.mean(), 'sz':0, 'st':0, 'sv':0}

# simulate complex task:
dfs = []
for i in range(n_subjects):
    df0 = simulate_data_GnG(zGo=params0['par_zGo'], zNoGoGo=params0['par_zNoGoGo'], a=params0['a'], vGo=params0['vGo'], vNoGoGo = params0['vNoGoGo'], 
                        t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'], 
                        condition=params0['cond'])
    df1 = simulate_data_GnG(zGo=params0['par_zGo'], zNoGoGo=params0['par_zNoGoGo'], a=params1['a'], vGo=params1['vGo'], vNoGoGo = params1['vNoGoGo'],
                        t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'],
                        condition=params1['cond'])
    df2 = simulate_data_GnG(zGo=params0['par_zGo'], zNoGoGo=params0['par_zNoGoGo'], a=params2['a'], vGo=params2['vGo'], vNoGoGo = params2['vNoGoGo'], 
                        t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'],
                        condition=params2['cond'])
    df = pd.concat((df0, df1, df2))
    df['subj_idx'] = i
    dfs.append(df)

# combine in one dataframe:
df_emp = pd.concat(dfs)
df_emp.loc[df_emp["response"]==0, 'rt'] = -1

seaborn_data_plotting(df_emp)

model_driftbias.plot_posteriors(['t', 'z_Intercept', 'z_C(condition)[T.NoGo - Go]', 'a_C(StimC)[T.20Hz]', 'a_Intercept'])