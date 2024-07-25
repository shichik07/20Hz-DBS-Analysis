# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:42:16 2023

@author: doex9445
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import hddm
from joblib import Parallel, delayed
from patsy import dmatrix
import random
import math


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


def simulate_data_GnG(a, zz, v, t, z, z_Go, vc, sv=0, sz=0, st=0, condition = 0, nr_trials1=100, nr_trials2=50, nr_trials3=50):

    """
    Simulates stim-coded data.
    """

    parameters1 = {'a':a + zz, 'v':v, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters2 = {'a':a + zz, 'v':v, 't':t, 'z':z_Go, 'sv':sv, 'sz': sz, 'st': st} # we set the parameter to 0.5 here because we expect go and nogo trials here to be equally likely
    parameters3 = {'a':a + zz, 'v':vc, 't':t, 'z':1-z, 'sv':sv, 'sz': sz, 'st': st}
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

def simulate_data_GnG_output(a, vNoGoGo, vStop, vGo, t, z, z_Go, sv=0, sz=0, st=0, condition = 0, nr_trials1=100, nr_trials2=50, nr_trials3=50):

    """
    Simulates stim-coded data.
    """

    parameters1 = {'a':a, 'v': vGo, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters2 = {'a':a, 'v': vNoGoGo, 't':t, 'z':z_Go, 'sv':sv, 'sz': sz, 'st': st} # we set the parameter to 0.5 here because we expect go and nogo trials here to be equally likely
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

def simulate_data_GnG_simple(a, zz, v, t, z, dc, sv, sz, st, condition, nr_trials2, nr_trials3):

    """
    Simulates stim-coded data.
    """

    parameters2 = {'a':a + zz, 'v':v+dc, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters3 = {'a':a + zz, 'v':v-dc, 't':t, 'z':1-z, 'sv':sv, 'sz': sz, 'st': st}
    
    df_sim2, params_sim2 = hddm.generate.gen_rand_data(params=parameters2, size=nr_trials2, subjs=1, subj_noise=0)
    df_sim2['condition'] = "Go"
    df_sim3, params_sim3 = hddm.generate.gen_rand_data(params=parameters3, size=nr_trials3, subjs=1, subj_noise=0)
    df_sim3['condition'] = "NoGo"
    df_sim = pd.concat((df_sim2, df_sim3))
    df_sim['correct'] = df_sim.apply(get_choice, 1)
    #df_sim['correct'] = df_sim['response'].astype(int)
    #df_sim['response'] = df_sim['bias_response'].astype(int)
    #df_sim['stimulus'] = np.array((np.array(df_sim['response']==1) & np.array(df_sim['correct']==1)) + (np.array(df_sim['response']==0) & np.array(df_sim['correct']==0)), dtype=int)
    #df_sim['StimC'] = condition
    #df_sim = df_sim.drop(columns=['bias_response'])

    return df_sim

def simulate_data_GnG_simple2(a, zz, v, t, z, dc, sv, sz, st,nr_trials, nr_subjects):

    """
    Simulates stim-coded data.
    """

    
    # session 2
    level1b = {'v':v, 'a':a, 't':t,'sv': 0, 'z':z, 'sz': 0, 'st': 0}
    level2b = {'v':v, 'a':a, 't':t,'sv': 0, 'z':1-z, 'sz': 0, 'st': 0}
    
    # simulate data
    data_b, params_b = hddm.generate.gen_rand_data({'go': level1b,
                                                    'nogo': level2b,},
                                                    size=nr_trials,
                                                    subjs=nr_subjects)

    return data_b


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
    box_plot = sns.boxplot(data = data[data['condition'] != 'NoGo'],
                           hue ='condition', y = 'rt', x = 'StimC', 
                           order=hue_ord, ax = ax1, palette= pal)
    #sns.despine(offset = 10, trim = True) # make the axis more exiting
    # median RT added to plot
    median_RT = data[data['condition'] != 'NoGo'].groupby(['StimC', 'condition'])['rt'].median().round(3)
    vertical_offset = data[data['condition'] != 'NoGo'].rt.median()* 0.05 # offset from median for display
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
    
def seaborn_data_plotting_simple(data):
    sns.set_style("whitegrid")
    pal = sns.color_palette('Set2')
    #set figure size
    sns.set(rc = {'figure.figsize': (10, 6)})
    #data = pd.DataFrame(df_emp_simple)
    data.condition.astype('category')
    #reset index
    data = data.reset_index()
    data["error"] = 1 - data["correct"]
    # set up subplots
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    # plot RT
    box_plot = sns.boxplot(data = data[data['condition'] != 'NoGo'],
                           hue ='condition', y = 'rt', 
                           ax = ax1, palette= pal)
        
    # plot accuracy values
    ax2 = fig.add_subplot(122)
    bar_plot = sns.barplot(data = data, x = 'condition', y = 'error',  
            ax = ax2, palette= pal)
    
    plt.close(2)
    plt.close(3)
    plt.tight_layout()
        
    
    

random.seed(852)
# settings
go_nogo = True # should we put all RTs for one choice alternative to NaN (go-no data)?
n_subjects = 16 #16
GoTrials = 150 #75
NoGoTrials =150 #25

#parameters1 = {'a':a, 'v':v-dc, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
#parameters2 = {'a':a + zz, 'v':v-dc, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
#parameters3 = {'a':a + zz, 'v':v+dc, 't':t, 'z':1-z, 'sv':sv, 'sz': sz, 'st': st}

# parameters:
params0 = {'cond':'OFF', 'v':1.3, 'a':2.0, 'zz':0,'t':0.3, 'z':0.6, 'z_Go': 0.5,'vc':-.5, 'sz':0, 'st':0, 'sv':0}
params1 = {'cond':'130Hz', 'v':1.3, 'a':2.0, 'zz':0,'t':0.3, 'z':0.6, 'z_Go': 0.5, 'vc':-.5, 'sz':0, 'st':0, 'sv':0}
params2 = {'cond':'20Hz', 'v':1.3, 'a':2.0, 'zz':0.8,'t':0.3, 'z':0.6, 'z_Go': 0.5, 'vc':-.5, 'sz':0, 'st':0, 'sv':0}

# simulate complex task:
dfs = []
for i in range(n_subjects):
    df0 = simulate_data_GnG(z=params0['z'], z_Go=params0['z_Go'], a=params0['a'], v=params0['v'], vc=params0['vc'],
                        t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'],
                        condition=params0['cond'], zz = params0['zz'])
    df1 = simulate_data_GnG(z=params1['z'], z_Go=params0['z_Go'], a=params1['a'], v=params1['v'], vc=params1['vc'],
                        t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'],
                        condition=params1['cond'], zz = params1['zz'])
    df2 = simulate_data_GnG(z=params2['z'], z_Go=params0['z_Go'], a=params2['a'], v=params2['v'], vc=params2['vc'],
                        t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'],
                        condition=params2['cond'], zz = params2['zz'])
    df = pd.concat((df0, df1, df2))
    df['subj_idx'] = i
    dfs.append(df)

# combine in one dataframe:
df_emp = pd.concat(dfs)
if go_nogo:
    df_emp.loc[df_emp["response"]==0, 'rt'] = -1
    
# Plot the accuracy data for Go and NoGo responses
seaborn_data_plotting(df_emp)

# next let us fit the model, see if we can recover the parameters

# Let us write a to do List: we need three Link functions:
    # 1) We need a Link function to implement the starting point bias for Go/NoGo distinction
    # 2) We need a Link function to implement the drift rate bias for the Go/NoGo distinction
    # 3) Worst of all - we need a link function to implement the difference in the decision threshold by Go/NoGo decision and Stimulation condition


# get an index for stop or go 
df_emp['stimulus'] = df_emp["condition"].apply(lambda x: 0 if x == "NoGo" else 1) # column to code Stop and Go trials
# get a sum contrast for the Go -NoGo vs Go condition 
df_emp['Contrast_Go'] = df_emp["condition"].apply(lambda x: 1 if x == "Go" else(-1 if x == "Go - NoGo" else 0))
# Then we need another column that shows us the NoGo trial condition (including the Go trials)
df_emp['ChangeCondition'] =df_emp["condition"].apply(lambda x: 1 if "NoGo" in x else 0)
# get a sum contrast for the NoGo vs Go condition 
df_emp['Contrast_GovsStop'] = df_emp["condition"].apply(lambda x: -1 if x == "NoGo" else 1)
 
# reset index, otherwise the link functions are a MESS
df_emp = df_emp.reset_index(drop = True)
mydata = df_emp

x = df_emp['Contrast_Go']

# first the link function for the starting point bias - which we took from the hddm website
def z_link_func(x, data=mydata):
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
#x = df_emp.stimulus
def v_link_func(x, data=mydata):
    stim = (np.asarray(dmatrix('0 + C(s, [[-1], [1]])',
                              {'s': data.stimulus.loc[x.index]},return_type='dataframe'))
    )
    return np.multiply(x.to_frame(), stim)

# def z_link_func(x, data=mydata):
#     stim = (np.asarray(dmatrix('0 + C(s, [[1], [-1]])',
#                                {'s': data.stimulus.loc[x.index]}))
#     )
#     return 1 / (1 + np.exp(-(x * stim)))

# Okay so it might be possible that we do not need a separate link function for the decision threshold
# Instead, what we are going to do is create a new variable that contains Go trials (Go) versus NoGo (NoGo + Go) trials
# and estimate the threshold of the interaction with the normal stimulus function. I do not think the threshold for Go versus NoGo should vary
# Only by the condition. So lets set this up and see if we can recover our simulated parameter
v_reg = {'model': 'v ~ 1 + Contrast_GovsStop', 'link_func': lambda v: v}
#v_reg = {'model': 'v ~ 1 + Contrast_GovsStop', 'link_func': lambda v: v_link_func}
z_reg = {'model': 'z ~ 1 + Contrast_Go', 'link_func': z_link_func}
#z_reg = {'model': 'z ~ 0 + C(condition)', 'link_func':  lambda z: z}
a_reg = {'model': "a ~ 1 + C(StimC)", 'link_func': lambda a: a}

# from patsy import dmatrices
# formula = "rt ~ 1 + C(StimC, Treatment('OFF')):C(change_con, Treatment(0))"
# y, X = dmatrices(formula, df_emp_sub, return_type='dataframe')
# X['StimC'] = df_emp_sub.StimC.copy()
# X['change_con'] = df_emp_sub.change_con.copy()


reg_model =  [z_reg, v_reg, a_reg]
sim_mod = hddm.HDDMRegressor(df_emp, 
                             reg_model, 
                             include=['z' , 'v', 'a'])

# find a good starting value
sim_mod.find_starting_values()
#get samples and discard a couple as burn ins
sim_mod.sample(600, burn = 100, dbname='GoNoGo_threshold_byc_Sim',db='pickle')
# save model
model = hddm.load('mymodel')
#check stats
staty = sim_mod.print_stats()

t, z, z_GoNoGo, v_I, v_C, a_130, a_20, a_OFF = sim_mod.nodes_db.loc[["t", "z_Intercept", "z_Contrast_Go", "v_Intercept", "v_Contrast_GovsStop",
                                                        "a_Intercept", "a_C(StimC)[T.20Hz]", "a_C(StimC)[T.OFF]"], 
                                                                  'node']


# use the mean as the parameter of interest
par_t = t.trace()
par_z = z.trace() - z_GoNoGo.trace() 
par_z_GoNoGo = z.trace() + z_GoNoGo.trace()
par_vGo = v_I.trace() + v_C.trace()
par_vStop = v_I.trace() - v_C.trace()
par_a130 = a_130.trace()
par_a20 = a_20.trace() + par_a130
par_aOFF =  a_OFF.trace() + par_a130


# next we check if the model parameter fit the data

# parameters:
n_subjects = 16
params0 = {'cond':'OFF', 'vGo': par_vGo.mean(), 'vNoGoGo':  par_vGo.mean(), 'vStop': par_vStop.mean(), 
           'a':par_aOFF.mean(), 't':par_t.mean(), 'z':par_z.mean(), 'z_Go': par_z_GoNoGo.mean(), 'sz':0, 'st':0, 'sv':0}
params1 = {'cond':'130Hz', 'vGo': par_vGo.mean(), 'vNoGoGo':  par_vGo.mean(), 'vStop': par_vStop.mean(), 
           'a':par_a130.mean(), 't':par_t.mean(), 'z':par_z.mean(), 'z_Go': par_z_GoNoGo.mean(), 'sz':0, 'st':0, 'sv':0}
params2 = {'cond':'20Hz', 'vGo': par_vGo.mean(), 'vNoGoGo':  par_vGo.mean(), 'vStop': par_vStop.mean(), 
           'a':par_a20.mean(), 't':par_t.mean(), 'z':par_z.mean(), 'z_Go': par_z_GoNoGo.mean(), 'sz':0, 'st':0, 'sv':0}

# simulate complex task:
dfs = []
for i in range(n_subjects):
    df0 = simulate_data_GnG_output(z=params0['z'], z_Go=params0['z_Go'], a=params0['a'], vGo=params0['vGo'], vNoGoGo = params0['vNoGoGo'], vStop = params0['vStop'],
                        t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'], 
                        condition=params0['cond'])
    df1 = simulate_data_GnG_output(z=params1['z'], z_Go=params0['z_Go'], a=params1['a'], vGo=params1['vGo'], vNoGoGo = params1['vNoGoGo'], vStop = params1['vStop'],
                        t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'],
                        condition=params1['cond'])
    df2 = simulate_data_GnG_output(z=params2['z'], z_Go=params0['z_Go'], a=params2['a'], vGo=params2['vGo'], vNoGoGo = params2['vNoGoGo'], vStop = params2['vStop'],
                        t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'],
                        condition=params2['cond'])
    df = pd.concat((df0, df1, df2))
    df['subj_idx'] = i
    dfs.append(df)

# combine in one dataframe:
df_sim = pd.concat(dfs)
df_sim.loc[df_sim["response"]==0, 'rt'] = -1
    
# Plot the accuracy data for Go and NoGo responses
seaborn_data_plotting(df_sim)

# Intercept_a20Hz, NoGo_change, Hz130Go, OFFGo, Hz130NoGo, OFFNoGO = sim_mod.nodes_db.loc[["a_Intercept",
#                                             "a_C(change_con, Treatment(0))[T.1]", 
#                                   "a_C(StimC, Treatment(1))[T.130Hz]:C(change_con, Treatment(0))[0]", 
#                                   "a_C(StimC, Treatment(1))[T.OFF]:C(change_con, Treatment(0))[0]",
#                                   "a_C(StimC, Treatment(1))[T.130Hz]:C(change_con, Treatment(0))[1]",
#                                   "a_C(StimC, Treatment(1))[T.OFF]:C(change_con, Treatment(0))[1]"
#                                   ], 
#                                            'node']

# # Next we construct the threshold off these parameters for each condition by stimulation
# hz20NoGo = Intercept_a20Hz.trace() + NoGo_change.trace()
# hz20Go = Intercept_a20Hz.trace() 
# hz130NoGo = Intercept_a20Hz.trace() + NoGo_change.trace() + Hz130NoGo.trace()
# hz130Go = Intercept_a20Hz.trace() + Hz130Go.trace()
# OffNoGo = Intercept_a20Hz.trace() + NoGo_change.trace() + OFFNoGO.trace()
# Off20Go = Intercept_a20Hz.trace() + OFFGo.trace()

                                
# # Plot and Store Distributions
# distributions = [hz20NoGo, hz20Go, hz130NoGo,
#                  hz130Go, OffNoGo,  Off20Go]

# # Plot densities using seaborn
# sns.set(style="whitegrid")
# for distribution in distributions:
#     sns.kdeplot(distribution, shade=True)

# plt.xlabel('decision threshold')
# plt.ylabel('Posterior probability')
# plt.title('Group mean posteriors of within-subject threshold effect Cued')
# plt.legend(['20Hz Change', '20Hz Go', '130Hz Change', "130Hz Go", "Off change", "Off Go"], 
#            loc='center left', bbox_to_anchor=(1, 0.5))
# plt.savefig('hddm_threshold_GoNoGO_Cued.pdf')

# ## This part needs to be rewritten

# print("P(20Hz > 130Hz Cued) = ", ((Hz20.trace() < Hz130.trace())).mean())
# print("P(20Hz > OFF Cued) = ", ((Hz20.trace() <  OFF.trace())).mean())
# print("P(OFF > 130Hz Cued) = ", ((OFF.trace() <  Hz130.trace())).mean())


# Hz1301, Hz201, OFF1 = md_srt_stim.nodes_db.loc[["a_C(Stim_verb, Treatment('OFF'))[130Hz]:C(change_con, Treatment(0))[0]", 
#                                   "a_C(Stim_verb, Treatment('OFF'))[20Hz]:C(change_con, Treatment(0))[0]", 
#                                   "a_C(Stim_verb, Treatment('OFF'))[OFF]:C(change_con, Treatment(0))[0]"], 'node']
# hddm.analyze.plot_posterior_nodes([Hz1301, Hz201, OFF1])
# plt.xlabel('decision threshold')
# plt.ylabel('Posterior probability')
# plt.title('Group mean posteriors of within-subject threshold effect Not Cued')
# plt.legend(['130Hz', '20Hz', 'OFF'], 
#            loc='center left', bbox_to_anchor=(1, 0.5))
# plt.savefig('hddm_threshold_GoNoGO_NoNCued.pdf')

# print("P(20Hz > 130Hz Non-Cued) = ", ((Hz201.trace() < Hz1301.trace())).mean())
# print("P(20Hz > OFF Non-Cued) = ", ((Hz201.trace() <  OFF1.trace())).mean())
# print("P(OFF > 130Hz Non-Cued) = ", ((OFF1.trace() <  Hz1301.trace())).mean())
