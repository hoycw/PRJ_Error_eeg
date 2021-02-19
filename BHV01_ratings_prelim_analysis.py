
# coding: utf-8

# In[1]:

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


# In[43]:

SBJ = sys.argv[1]#raw_input('Enter SBJ ID to process:')#'EEG15'


# In[44]:

prj_dir = '/Volumes/hoycw_clust/PRJ_Error_eeg/'#'/Users/sheilasteiner/Desktop/Knight_Lab/PRJ_Error_eeg/'
results_dir = prj_dir+'results/'
fig_type = '.png'
data_dir = prj_dir+'data/'
sbj_dir  = data_dir+SBJ+'/'


# ### Load paradigm parameters

# In[45]:

prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars.pkl')
with open(prdm_fname, 'rb') as f:
    prdm = pickle.load(f)


# ### Load Log Info

# In[46]:

behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav.csv')
data = pd.read_csv(behav_fname)


# In[47]:

# Remove second set of training trials in restarted runs (EEG12, EEG24, EEG25)
if len(data[(data['Trial']==0) & (data['Block']==-1)])>1:
    train_start_ix = data[(data['Trial']==0) & (data['Block']==-1)].index
    train_ix = [ix for ix in data.index if data.loc[ix,'Block']==-1]
    later_ix = [ix for ix in data.index if ix >= train_start_ix[1]]
    data = data.drop(set(later_ix).intersection(train_ix))
    data = data.reset_index()


# In[48]:

# Change block numbers on EEG12 to not overlap
if SBJ=='EEG12':
    b4_start_ix = data[(data['Trial']==0) & (data['Block']==4)].index
    for ix in range(b4_start_ix[1]):
        if data.loc[ix,'Block']!=-1:
            data.loc[ix,'Block'] = data.loc[ix,'Block']-4


# In[49]:

# Label post-correct (PC), post-error (PE) trials
data['PE'] = [False for _ in range(len(data))]
for ix in range(len(data)):
    # Exclude training data and first trial of the block
    if (data.loc[ix,'Block']!=-1) and (data.loc[ix,'Trial']!=0):
        if data.loc[ix-1,'Hit']==0:
            data.loc[ix,'PE'] = True


# In[50]:

# pd.set_option('max_rows', 75)
# data[data['Block']==3]


# # Add specific analysis computations

# In[51]:

# Find middle of blocks to plot accuracy
block_start_ix = data[data['Trial']==0].index
if SBJ=='EP11':#deal with missing BT_T0
    block_mid_ix = [ix+prdm['n_trials']/2 for ix in block_start_ix]
else:
    block_mid_ix = [ix+prdm['n_trials']/2 for ix in block_start_ix[1:]]

# Add in full_vis + E/H training: 0:4 + 5:19 = 10; 20:34 = 27.5 
block_mid_ix.insert(0,np.mean([prdm['n_examples']+prdm['n_training'],
         prdm['n_examples']+2*prdm['n_training']]))   #examples
block_mid_ix.insert(0,np.mean([0, prdm['n_examples']+prdm['n_training']]))
#easy training (would be 12.5 if splitting examples/train)


# In[52]:

# Compute accuracy per block
accuracy = data['Hit'].groupby([data['Block'],data['Condition']]).mean()
acc_ITI = data['Hit'].groupby([data['ITI type'],data['Condition']]).mean()
for ix in range(len(data)):
    data.loc[ix,'Accuracy'] = accuracy[data.loc[ix,'Block'],data.loc[ix,'Condition']]
    data.loc[ix,'Acc_ITI'] = acc_ITI[data.loc[ix,'ITI type'],data.loc[ix,'Condition']]


# In[53]:

# Break down by post-long and post-short trials
data['postlong'] = [False if ix==0 else True if data['RT'].iloc[ix-1]>1 else False for ix in range(len(data))]

# Compute change in RT
data['dRT'] = [0 for ix in range(len(data))]
for ix in range(len(data)-1):
    data.loc[ix+1,'dRT'] = data.loc[ix+1,'RT']-data.loc[ix,'RT']


# In[54]:

# Grab rating data to plot
rating_trial_idx = [True if rating != -1 else False for rating in data['Rating']]
rating_data = data['Rating'][rating_trial_idx]


# # Plot Full Behavior Across Dataset

# In[55]:

# Accuracy, Ratings, and Tolerance
f, ax1 = plt.subplots()
x = range(len(data))
plot_title = '{0} Tolerance and Accuracy: easy={1:0.3f}; hard={2:0.3f}'.format(
                SBJ, data[data['Condition']=='easy']['Hit'].mean(),
                data[data['Condition']=='hard']['Hit'].mean())
    
colors = {'easy': [0.5, 0.5, 0.5],#[c/255 for c in [77,175,74]],
          'hard': [1, 1, 1],#[c/255 for c in [228,26,28]],
          'accuracy': 'k'}#[c/255 for c in [55,126,184]]}
scat_colors = {'easy': [1,1,1],#[c/255 for c in [77,175,74]],
          'hard': [0,0,0]}
accuracy_colors = [scat_colors[accuracy.index[ix][1]] for ix in range(len(accuracy))]
#scale = {'Hit Total': np.max(data['Tolerance'])/np.max(data['Hit Total']),
#         'Score Total': np.max(data['Tolerance'])/np.max(data['Score Total'])}

# Plot Tolerance Over Time
ax1.plot(data['Tolerance'],'b',label='Tolerance')
ax1.plot(x,[prdm['tol_lim'][0] for _ in x],'b--')
ax1.plot(x,[prdm['tol_lim'][1] for _ in x],'b--')
ax1.set_ylabel('Target Tolerance (s)', color='b')
ax1.tick_params('y', colors='b')
ax1.set_xlim([0,len(data)])
ax1.set_ylim([0, 0.41])
ax1.set_facecolor('white')
ax1.grid(False)

# Plot Accuracy per Block
ax2 = ax1.twinx()
# ax2.plot(data['Hit Total']/np.max(data['Hit Total']),'k',label='Hit Total')
ax2.fill_between(x, 1, 0, where=data['Condition']=='easy',
                facecolor=colors['easy'], alpha=0.3)#, label='hard')
ax2.fill_between(x, 1, 0, where=data['Condition']=='hard',
                facecolor=colors['hard'], alpha=0.3)#, label='easy')
ax2.scatter(block_mid_ix, accuracy, s=50, c=accuracy_colors,
           edgecolors='k', linewidths=1)#colors['accuracy'])#,linewidths=2)
ax2.scatter(rating_data.index.values, rating_data.values/100, s=25, c=[1, 0, 0])
ax2.set_ylabel('Accuracy', color=colors['accuracy'])
ax2.tick_params('y', colors=colors['accuracy'])
ax2.set_xlabel('Trials')
ax2.set_xlim([0,len(data)])
ax2.set_ylim([0, 1])
ax2.set_facecolor('white')
ax2.grid(False)

plt.title(plot_title)

plt.savefig(results_dir+'BHV/ratings_tolerance/'+SBJ+'_tolerance'+fig_type)
plt.close()


# In[ ]:


