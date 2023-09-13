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

def simulate_data(a, v, t, z, dc, sv=0, sz=0, st=0, condition=0, nr_trials1=100, nr_trials2=50):

    """
    Simulates stim-coded data.
    """

    parameters1 = {'a':a, 'v':v+dc, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters2 = {'a':a, 'v':v-dc, 't':t, 'z':1-z, 'sv':sv, 'sz': sz, 'st': st}
    df_sim1, params_sim1 = hddm.generate.gen_rand_data(params=parameters1, size=nr_trials1, subjs=1, subj_noise=0)
    df_sim1['condition'] = 'present'
    df_sim2, params_sim2 = hddm.generate.gen_rand_data(params=parameters2, size=nr_trials2, subjs=1, subj_noise=0)
    df_sim2['condition'] = 'absent'
    df_sim = pd.concat((df_sim1, df_sim2))
    df_sim['bias_response'] = df_sim.apply(get_choice, 1)
    df_sim['correct'] = df_sim['response'].astype(int)
    df_sim['response'] = df_sim['bias_response'].astype(int)
    df_sim['stimulus'] = np.array((np.array(df_sim['response']==1) & np.array(df_sim['correct']==1)) + (np.array(df_sim['response']==0) & np.array(df_sim['correct']==0)), dtype=int)
    df_sim['condition'] = condition
    df_sim = df_sim.drop(columns=['bias_response'])

    return df_sim

def simulate_data_GnG(a, zz, v, t, z, dc, sv=0, sz=0, st=0, condition = 0, nr_trials1=50, nr_trials2=25, nr_trials3=25):

    """
    Simulates stim-coded data.
    """

    parameters1 = {'a':a, 'v':v+dc, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters2 = {'a':a + zz, 'v':v+dc, 't':t, 'z':z, 'sv':sv, 'sz': sz, 'st': st}
    parameters3 = {'a':a + zz, 'v':v-dc, 't':t, 'z':1-z, 'sv':sv, 'sz': sz, 'st': st}
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

def fit_subject(data, quantiles):

    """
    Simulates stim-coded data.
    """

    subj_idx = np.unique(data['subj_idx'])
    m = hddm.HDDMStimCoding(data, stim_col='stimulus', split_param='v', drift_criterion=True, bias=True, p_outlier=0,
                            depends_on={'v':'condition', 'a':'condition', 't':'condition', 'z':'condition', 'dc':'condition', })
    m.optimize('gsquare', quantiles=quantiles, n_runs=8)
    res = pd.concat((pd.DataFrame([m.values], index=[subj_idx]), pd.DataFrame([m.bic_info], index=[subj_idx])), axis=1)
    return res

def summary_plot(df_group, df_sim_group=None, quantiles=[0, 0.1, 0.3, 0.5, 0.7, 0.9,], xlim=None):

    """
    Generates a
    """

    nr_subjects = len(np.unique(df_group['subj_idx']))

    fig = plt.figure(figsize=(10,nr_subjects*2))
    plt_nr = 1
    for s in np.unique(df_group['subj_idx']):
        df = df_group.copy().loc[(df_group['subj_idx']==s),:]
        df_sim = df_sim_group.copy().loc[(df_sim_group['subj_idx']==s),:]
        df['rt_acc'] = df['rt'].copy()
        df.loc[df['correct']==0, 'rt_acc'] = df.loc[df['correct']==0, 'rt_acc'] * -1
        df['rt_resp'] = df['rt'].copy()
        df.loc[df['response']==0, 'rt_resp'] = df.loc[df['response']==0, 'rt_resp'] * -1
        df_sim['rt_acc'] = df_sim['rt'].copy()
        df_sim.loc[df_sim['correct']==0, 'rt_acc'] = df_sim.loc[df_sim['correct']==0, 'rt_acc'] * -1
        df_sim['rt_resp'] = df_sim['rt'].copy()
        df_sim.loc[df_sim['response']==0, 'rt_resp'] = df_sim.loc[df_sim['response']==0, 'rt_resp'] * -1
        max_rt = np.percentile(df_sim.loc[~np.isnan(df_sim['rt']), 'rt'], 99)
        bins = np.linspace(-max_rt,max_rt,30)
        # rt distributions correct vs error:
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
        N, bins, patches = ax.hist(df.loc[:, 'rt_acc'], bins=bins,
                                   density=True, color='green', alpha=0.5)


        for bin_size, bin, patch in zip(N, bins, patches):
            if bin < 0:
                plt.setp(patch, 'facecolor', 'r')
        if df_sim is not None:
            ax.hist(df_sim.loc[:, 'rt_acc'], bins=bins, density=True,
                    histtype='step', color='k', alpha=1, label=None)
        ax.set_title('P(correct)={}'.format(round(df.loc[:, 'correct'].mean(), 3),))
        ax.set_xlabel('RT (s)')
        ax.set_ylabel('Trials (prob. dens.)')
        plt_nr += 1

        # condition accuracy plots:
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
        df.loc[:,'rt_bin'] = pd.qcut(df['rt'], quantiles, labels=False)
        d = df.groupby(['rt_bin']).mean().reset_index()
        ax.errorbar(d.loc[:, "rt"], d.loc[:, "correct"], fmt='-o', color='orange', markersize=10)
        if df_sim is not None:
            df_sim.loc[:,'rt_bin'] = pd.qcut(df_sim['rt'], quantiles, labels=False)
            d = df_sim.groupby(['rt_bin']).mean().reset_index()
            ax.errorbar(d.loc[:, "rt"], d.loc[:, "correct"], fmt='x', color='k', markersize=6)
        if xlim:
            ax.set_xlim(xlim)
        ax.set_ylim(0, 1.25)
        ax.set_title('Conditional accuracy')
        ax.set_xlabel('RT (quantiles)')
        ax.set_ylabel('P(correct)')
        plt_nr += 1

        # rt distributions response 1 vs 0:
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
        if np.isnan(df['rt']).sum() > 0:
            # some initial computations
            bar_width = 1
            fraction_yes = df['response'].mean()
            fraction_yes_sim = df_sim['response'].mean()
            no_height = (1 - fraction_yes) / bar_width
            no_height_sim = (1 - fraction_yes_sim) / bar_width

            hist, edges = np.histogram(df.loc[:, 'rt_resp'], bins=bins, density=True,)
            hist = hist * fraction_yes
            hist_sim, edges_sim = np.histogram(df_sim.loc[:, 'rt_resp'], bins=bins, density=True,)
            hist_sim = hist_sim * fraction_yes_sim

            # Add histogram from go choices
            # ground truth
            ax.bar(edges[:-1], hist, width=np.diff(edges)[0], align='edge',
                   color='magenta', alpha=0.5, linewidth=0,)
            # simulations
            ax.step(edges_sim[:-1] + np.diff(edges)[0], hist_sim, color='black', lw=1)

            # Add bar for the no-go choices (on the negative rt scale)
            # This just illustrates the probability of no-go choices

            # ground truth
            ax.bar(x=-1.5, height=no_height, width=bar_width, alpha=0.5, color='cyan', align='center')

            # simulations
            ax.hlines(y=no_height_sim, xmin=-2, xmax=-1, lw=0.5, colors='black',)
            ax.vlines(x=-2, ymin=0, ymax=no_height_sim, lw=0.5, colors='black')
            ax.vlines(x=-1, ymin=0, ymax=no_height_sim, lw=0.5, colors='black')
        else:
            N, bins, patches = ax.hist(df.loc[:, 'rt_resp'], bins=bins,
                                   density=True, color='magenta', alpha=0.5)
            for bin_size, bin, patch in zip(N, bins, patches):
                if bin < 0:
                    plt.setp(patch, 'facecolor', 'cyan')
            ax.hist(df_sim.loc[:, 'rt_resp'], bins=bins, density=True,
                    histtype='step', color='k', alpha=1, label=None)

        ax.set_title('P(bias)={}'.format(round(df.loc[:, 'response'].mean(), 3),))
        ax.set_xlabel('RT (s)')
        ax.set_ylabel('Trials (prob. dens.)')
        plt_nr += 1

        # condition response plots:
        ax = fig.add_subplot(nr_subjects,4,plt_nr)
        df.loc[:,'rt_bin'] = pd.qcut(df['rt'], quantiles, labels=False)
        d = df.groupby(['rt_bin']).mean().reset_index()
        ax.errorbar(d.loc[:, "rt"], d.loc[:, "response"], fmt='-o', color='orange', markersize=10)
        if df_sim is not None:
            df_sim.loc[:,'rt_bin'] = pd.qcut(df_sim['rt'], quantiles, labels=False)
            d = df_sim.groupby(['rt_bin']).mean().reset_index()
            ax.errorbar(d.loc[:, "rt"], d.loc[:, "response"], fmt='x', color='k', markersize=6)
        if xlim:
            ax.set_xlim(xlim)
        ax.set_ylim(0,1.25)
        ax.set_title('Conditional response')
        ax.set_xlabel('RT (quantiles)')
        ax.set_ylabel('P(bias)')
        plt_nr += 1

    sns.despine(offset=3, trim=True)
    plt.tight_layout()

    return fig

random.seed(852)
# settings
go_nogo = True # should we put all RTs for one choice alternative to NaN (go-no data)?
n_subjects = 16
GoTrials = 75
NoGoTrials =25

# parameters:
params0 = {'cond':'OFF', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':0.7, 'dc':-0.2, 'sz':0, 'st':0, 'sv':0}
params1 = {'cond':'130Hz', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':0.7, 'dc':-0.2, 'sz':0, 'st':0, 'sv':0}
params2 = {'cond':'20Hz', 'v':.5, 'a':1.5, 'zz':0.5,'t':0.3, 'z':.7, 'dc':-0.2, 'sz':0, 'st':0, 'sv':0}


# # parameters:
# params0 = {'cond':'OFF', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':0.7, 'dc':0, 'sz':0, 'st':0, 'sv':0}
# params1 = {'cond':'130Hz', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':0.7, 'dc':0, 'sz':0, 'st':0, 'sv':0}
# params2 = {'cond':'20Hz', 'v':.5, 'a':1.5, 'zz':0,'t':0.3, 'z':.7, 'dc':0, 'sz':0, 'st':0, 'sv':0}


# simulate:
dfs = []
for i in range(n_subjects):
    df0 = simulate_data_GnG(z=params0['z'], a=params0['a'], v=params0['v'], dc=params0['dc'],
                        t=params0['t'], sv=params0['sv'], st=params0['st'], sz=params0['sz'],
                        condition=params0['cond'], zz = params0['zz'])
    df1 = simulate_data_GnG(z=params1['z'], a=params1['a'], v=params1['v'], dc=params1['dc'],
                        t=params1['t'], sv=params1['sv'], st=params1['st'], sz=params1['sz'],
                        condition=params1['cond'], zz = params1['zz'])
    df2 = simulate_data_GnG(z=params2['z'], a=params2['a'], v=params2['v'], dc=params2['dc'],
                        t=params2['t'], sv=params2['sv'], st=params2['st'], sz=params2['sz'],
                        condition=params2['cond'], zz = params2['zz'])
    df = pd.concat((df0, df1, df2))
    df['subj_idx'] = i
    dfs.append(df)

# combine in one dataframe:
df_emp = pd.concat(dfs)
if go_nogo:
    df_emp.loc[df_emp["response"]==0, 'rt'] = 0
    
# Plot the accuracy data for Go and NoGo responses

df = df_emp.copy()
# only plot correctly answered trials
fig = plt.figure(figsize=(20, 20))
# rt distributions correct vs error: df = df_group.copy().loc[(df_group['subj_idx']==s),:]
df['rt_acc'] = df['rt'].copy()
df.loc[df['correct']==0, 'rt_acc'] = df.loc[df['correct']==0, 'rt_acc'] * -1
df['rt_resp'] = df['rt'].copy()
df.loc[df['response']==0, 'rt_resp'] = df.loc[df['response']==0, 'rt_resp'] * -1
max_rt = np.percentile(df.loc[~np.isnan(df['rt']), 'rt'], 99)
bins = np.linspace(-max_rt,max_rt,30)
# rt distributions correct vs error:
ax = fig.add_subplot(5,6,1)
N, bins, patches = ax.hist(df.loc[:, 'rt_acc'], bins=bins,
                   density=True, color='green', alpha=0.5)


for bin_size, bin, patch in zip(N, bins, patches):
     if bin < 0:
         plt.setp(patch, 'facecolor', 'r')

ax.set_title('P(correct)={}'.format(round(df.loc[:, 'correct'].mean(), 3),))
ax.set_xlabel('RT (s)')
ax.set_ylabel('Trials (prob. dens.)')

# condition accuracy plots:
quantiles=[0, 0.1, 0.3, 0.5, 0.7, 0.9,]
ax = fig.add_subplot(5,5,2)
df.loc[:,'rt_bin'] = pd.qcut(df['rt'], quantiles, labels=False)
d = df.groupby(['rt_bin']).mean().reset_index()
ax.errorbar(d.loc[:, "rt"], d.loc[:, "correct"], fmt='-o', color='orange', markersize=10)
ax.set_ylim(0, 1.25)
ax.set_title('Conditional accuracy')
ax.set_xlabel('RT (quantiles)')
ax.set_ylabel('P(correct)')

# add accuracy by Go and NoGo responses
# ax = fig.add_subplot(5,5,3)
# d2 = df.copy()
# d2['Go'] =d2['condition'].apply(lambda x: 1 if 'NoGo' in x else 0)
# d2 = d2.groupby(['condition']).mean().reset_index()


colors = ['green', 'blue', 'red']
# bars = ax.bar(d2.loc[:, 'condition'], d2.loc[:, 'correct'], alpha=0.5)

# for bar, color in zip(bars, colors):
#     bar.set_color(color)

# ax.set_title('Accuracy by Stimulus')
# ax.set_xlabel('Response')
# ax.set_ylabel('P(correct)')

# add accuracy by Go and NoGo and NoGo-Go responses
ax = fig.add_subplot(5,5,3)
d3 = df.copy()
d3 = d3.groupby(['StimC','condition']).mean().reset_index()
sns.barplot(data = d3, x = "StimC", y = "correct", hue = "condition", palette = colors, alpha = 0.5 )
plt.legend([], [], frameon=False)
ax = fig.add_subplot(5,5,4)
sns.barplot(data = d3, x = "StimC", y = "rt", hue = "condition", palette = colors, alpha = 0.5 )
plt.legend([], [], frameon=False)

# next let us fit the model, see if we can recover the parameters

# Let us write a to do List: we need three Link functions:
    # 1) We need a Link function to implement the starting point bias for Go/NoGo distinction
    # 2) We need a Link function to implement the drift rate bias for the Go/NoGo distinction
    # 3) Worst of all - we need a link function to implement the difference in the decision threshold by Go/NoGo decision and Stimulation condition
# first create a new column with the Go versus NoGo condition
df_emp['change_con'] = df_emp["condition"].apply(lambda x: 1 if "NoGo" in x else 0)

# we want another condition for the interaction effect
# Define the strings to search for in each column
string1 = '20Hz'
string2 = 'NoGo'

# Define the function to apply on each row
def check_strings(row):
    if string1 in row['StimC'] and string2 in row['condition']:
        return 1
    else:
        return 0

# Apply the function to create a new column
df_emp['thresh'] = df_emp.apply(lambda row: check_strings(row), axis=1)

# okay we need to streamline the size of the dataframe. First we drop all columns that do not matter for the analysis
df_emp_sub = df_emp.copy().loc[:, ['rt', 'response', 'subj_idx', 'correct', 'stimulus', 'StimC', 'change_con', 'thresh']]
# change the stimC variable to integers
# Define the mapping dictionary for recoding
mapping = {
    'OFF': 1,
    '130Hz': 2,
    '20Hz': 3
}
 
# reset index, otherwise the link functions are a MESS
df_emp_sub = df_emp_sub.reset_index()

 
# first the link function for the starting point bias - which we took from the hddm website
def z_link_func(x, data=df_emp_sub):
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

def v_link_func(x, data=df_emp_sub):
    stim = (np.asarray(dmatrix('0 + C(s,[[1],[-1]])', {'s':data.stimulus.loc[x.index]})))
    return x * stim

# Okay so it might be possible that we do not need a separate link function for the decision threshold
# Instead, what we are going to do is create a new variable that contains Go trials (Go) versus NoGo (NoGo + Go) trials
# and estimate the threshold of the interaction with the normal stimulus function. I do not think the threshold for Go versus NoGo should vary
# Only by the condition. So lets set this up and see if we can recover our simulated parameter
v_reg = {'model': 'v ~ 1 + C(StimC)', 'link_func': v_link_func}
z_reg = {'model': 'z ~ 1 + C(StimC)', 'link_func': z_link_func}

a_reg = {'model': "a ~ 1 + C(StimC, Treatment(1)):C(change_con, Treatment(0))", 'link_func': lambda x: x}

# from patsy import dmatrices
# formula = "rt ~ 1 + C(StimC, Treatment('OFF')):C(change_con, Treatment(0))"
# y, X = dmatrices(formula, df_emp_sub, return_type='dataframe')
# X['StimC'] = df_emp_sub.StimC.copy()
# X['change_con'] = df_emp_sub.change_con.copy()


reg_model =  [z_reg, a_reg]
sim_mod = hddm.HDDMRegressor(df_emp_sub, 
                             reg_model, 
                             include=['v', 'a', 'z', 't'])

sim_mod = hddm.HDDMRegressor(df_emp_sub, 
                             {'model': "a ~ 1 + C(StimC, Treatment(1)):C(change_con, Treatment(0))", 'link_func': lambda x: x}, 
                             include=['v', 'a', 'z', 't'])
# find a good starting value
sim_mod.find_starting_values()
#get samples and discard a couple as burn ins
sim_mod.sample(2000, burn = 200, dbname='GoNoGo_threshold_byc_Sim',db='pickle')

staty = sim_mod.print_stats()

Intercept_a20Hz, NoGo_change, Hz130Go, OFFGo, Hz130NoGo, OFFNoGO = sim_mod.nodes_db.loc[["a_Intercept",
                                            "a_C(change_con, Treatment(0))[T.1]", 
                                  "a_C(StimC, Treatment(1))[T.130Hz]:C(change_con, Treatment(0))[0]", 
                                  "a_C(StimC, Treatment(1))[T.OFF]:C(change_con, Treatment(0))[0]",
                                  "a_C(StimC, Treatment(1))[T.130Hz]:C(change_con, Treatment(0))[1]",
                                  "a_C(StimC, Treatment(1))[T.OFF]:C(change_con, Treatment(0))[1]"
                                  ], 
                                           'node']

# Next we construct the threshold off these parameters for each condition by stimulation
hz20NoGo = Intercept_a20Hz.trace() + NoGo_change.trace()
hz20Go = Intercept_a20Hz.trace() 
hz130NoGo = Intercept_a20Hz.trace() + NoGo_change.trace() + Hz130NoGo.trace()
hz130Go = Intercept_a20Hz.trace() + Hz130Go.trace()
OffNoGo = Intercept_a20Hz.trace() + NoGo_change.trace() + OFFNoGO.trace()
Off20Go = Intercept_a20Hz.trace() + OFFGo.trace()

                                
# Plot and Store Distributions
distributions = [hz20NoGo, hz20Go, hz130NoGo,
                 hz130Go, OffNoGo,  Off20Go]

# Plot densities using seaborn
sns.set(style="whitegrid")
for distribution in distributions:
    sns.kdeplot(distribution, shade=True)

plt.xlabel('decision threshold')
plt.ylabel('Posterior probability')
plt.title('Group mean posteriors of within-subject threshold effect Cued')
plt.legend(['20Hz Change', '20Hz Go', '130Hz Change', "130Hz Go", "Off change", "Off Go"], 
           loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('hddm_threshold_GoNoGO_Cued.pdf')

## This part needs to be rewritten

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