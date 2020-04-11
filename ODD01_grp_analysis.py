from __future__ import division
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as io
import pickle

import scipy.stats

SBJ = sys.argv[1]

if os.path.exists('/Volumes/hoycw_clust/PRJ_Error_eeg/'):
    print('yes!')
    prj_dir='/Volumes/hoycw_clust/PRJ_Error_eeg/'
else:
    prj_dir = '/Users/sheilasteiner/Desktop/Knight_Lab/PRJ_Error_eeg/'
results_dir = prj_dir+'results/'
fig_type = '.png'
data_dir = prj_dir+'data/'
sbj_dir  = data_dir+SBJ+'/'

# ### Load paradigm parameters
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_odd_prdm_vars.pkl')
with open(prdm_fname, 'rb') as f:
    prdm = pickle.load(f)

#Load Behavior
behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav_oddball.csv')
data = pd.read_csv(behav_fname)

#Initialize Variables
block_range = np.arange(np.max(data['Block']+1))
condition_labels = ['tar', 'std', 'odd']
missed_trials = np.zeros((3,1))
accurate_ratio = np.zeros((3,3))

# Compute accuracy per condition
accurate_trials = np.zeros((3,1))
for cond_ix in range(len(condition_labels)):
    cond_trials = data[data['Condition'] == condition_labels[cond_ix]].index # find all the trials of a given condition
    for x in cond_trials:
        if condition_labels[cond_ix] == 'tar':
            if data.loc[x,'Hit'] == 0: # if got it right (Hit), increase tally
                missed_trials[cond_ix] +=1
        else:
            if data.loc[x,'Miss'] == 1: # count 1 - number of misses/total as measure of accuracy
                missed_trials[cond_ix] +=1

data_all = data

# Exclude: Training/Examples, first trial of each block
data = data[(data['Block']!=-1) & (data['ITI']>0)]
condition_titles = ['Target', 'Standard', 'Oddball']
return missed_trials
