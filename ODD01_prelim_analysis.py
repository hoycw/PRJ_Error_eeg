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
accurate_trials = np.zeros((3,1))
accurate_ratio = np.zeros((3,3))

# Compute accuracy per condition
for ix in block_range:
    accurate_trials = np.zeros((3,1))
    for cond_ix in range(len(condition_labels)):
        block_data = data[data['Block'] == ix] #get all the trials in a certain block
        cond_trials = block_data[block_data['Condition'] == condition_labels[cond_ix]].index # find all the trials of a given condition
        for x in cond_trials:
            if condition_labels[cond_ix] == 'tar':
                if data.loc[x,'Hit'] == 1: # if got it right (Hit), increase tally
                    accurate_trials[cond_ix] +=1
            else:
                if data.loc[x,'Miss'] == 0: # count 1 - number of misses/total as measure of accuracy
                    accurate_trials[cond_ix] +=1
        accurate_ratio[cond_ix,ix] = (accurate_trials[cond_ix]/np.size(cond_trials))# add the ratio of right/all (1 value for each block and each condition)
data_all = data

# Exclude: Training/Examples, first trial of each block
data = data[(data['Block']!=-1) & (data['ITI']>0)]
condition_titles = ['Target', 'Standard', 'Oddball']
# plot for each block the number correct, separated by condition
f, axes = plt.subplots(1,3)
for index in range(len(condition_titles)):
    sns.lineplot(block_range, accurate_ratio[index,:], ax=axes[index], markers = 'True', marker = "o")
    plt.subplots_adjust(top=0.8,wspace=0.8)
    axes[index].set_xticks([0,1,2])
    axes[index].set_xlabel('Block Number')
    axes[index].set_ylabel('Accuracy Rate')
    axes[index].set_ylim(0, 1.05)
    axes[index].set_title(condition_titles[index])

f.suptitle(SBJ + ' Condition and Accuracy in Oddball Task') # can also get the figure from plt.gcf()
if os.path.isdir(results_dir + 'BHV/ODD/accuracy/') == False:
    os.makedirs(results_dir + 'BHV/ODD/accuracy/')
plt.savefig(results_dir+'BHV/ODD/accuracy/'+SBJ+'_acc_condition'+fig_type)

