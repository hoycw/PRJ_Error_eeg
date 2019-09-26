
# coding: utf-8

# ### This version saves paradigm parameters and does not include some summary variables that can be computed from the main info:
#     reversal, post_err, score_total, hit_total, block_acc

# In[1]:

# get_ipython().magic(u'matplotlib inline')
import sys 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as io
import pickle


# In[2]:

SBJ = 'EEG01'


# In[3]:

#prj_dir = '/Volumes/hoycw_clust/PRJ_Error_eeg/'
prj_dir = '/Users/sheilasteiner/Desktop/Knight_Lab/PRJ_Error_eeg/'
results_dir = prj_dir+'results/'
fig_type = '.png'
data_dir = prj_dir+'data/'
sbj_dir  = data_dir+SBJ+'/'
# paths = {'Rana': '/Users/colinhoy/Code/PRJ_Error/data/logs/',
#          'Adi': '/Users/colinhoy/Code/PRJ_Error/data/logs/',
#          'IR57': '/Users/colinhoy/Code/PRJ_Error/data/logs/'}
logs = {'EEG01': 'eeg01_oddball_log_20190712145416.txt'
       }


# ### Load SBJ Data

# In[4]:

log_filename = os.path.join(sbj_dir,'00_raw',logs[SBJ])

log_file = open(log_filename,'r')
log = log_file.readlines()
log_file.close()


# ### Process Version-Specific Parameters

# In[5]:

prdm = {}
for line in log:
    # Script version
    if line.find('paradigm_name =')!=-1:
        prdm['prdm_name'] = line[line.find('= ')+2:line.find('\n')]
    if line.find('paradigm_version =')!=-1:
        prdm['prdm_version'] = line[line.find('= ')+2:line.find('\n')]
    if line.find('paradigm_type =')!=-1:
        prdm['prdm_type'] = line[line.find('= ')+2:line.find('\n')]

    # Design variables
    if line.find('use_rtbox =')!=-1:
        prdm['use_rtbox'] = bool(line[line.find('= ')+2:line.find('\n')])
    if line.find('n_blocks')!=-1:
        prdm['n_blocks'] = int(line[line.find('=')+2:])
    if line.find('n_trials')!=-1:
        prdm['n_trials'] = int(line[line.find('=')+2:])
    if line.find('n_training')!=-1:
        prdm['n_training'] = int(line[line.find('=')+2:])
    if line.find('block_order =')!=-1:
        prdm['block_order'] = [int(string)                     for string in line[line.find('[')+1:line.find(']')].split(' ')]
    if line.find('oddball order =')!=-1:
        prdm['odd_order'] = [int(string)                     for string in line[line.find('[')+1:line.find(']')].split(',')]
    if line.find('stim_dur =')!=-1:
        prdm['stim_dur'] = float(line[line.find('= ')+2:line.find('\n')])
    if line.find('visual_assignment =')!=-1:
        prdm['visual_assign'] = [int(string)                     for string in line[line.find('[')+1:line.find(']')].split(',')]
    
    # ITIs and boundaries between them
    if line.find('ITIs')!=-1:
        prdm['ITIs'] = [float(string)                     for string in line[line.find('[')+1:line.find(']')].split(',')]
        ITI_bounds = np.mean([prdm['ITIs'][:-1], prdm['ITIs'][1:]],0)
        


# ### Save paradigm parameters

# In[6]:

# Python readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_odd_prdm_vars.pkl')
with open(prdm_fname, 'wb') as f:
    pickle.dump(prdm, f, pickle.HIGHEST_PROTOCOL)
# MATLAB Readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_odd_prdm_vars.mat')
io.savemat(prdm_fname,prdm)


# In[7]:

print ('paradigm: ', prdm['prdm_name'], ' v', prdm['prdm_version'], ' type: ', prdm['prdm_type'])
print
print ('Stimulus duration: ', prdm['stim_dur'], 's')
print
print ('n_blocks: ', prdm['n_blocks'])
print ('n_trials/block: ', prdm['n_trials'])
print ('n_training: ', prdm['n_training'])
print
print ('block_order: ', prdm['block_order'])
print ('oddball_order: ', prdm['odd_order'])
# ITI_bounds = [np.mean(a,b) for a, b in zip(ITIs[:-1],ITIs[1:])]
print ('ITIs:',prdm['ITIs'], ITI_bounds)


# ### Extract Trial Info

# In[8]:

resp_lines = [line for line in log if line.find('Outcome =')!=-1]
data = pd.DataFrame({'Block': [line[line.find('B')+1] for line in resp_lines],
                     'Trial': [int(line[line.find('_T')+2:line.find(':')]) for line in resp_lines],
                     'Hit': [line.count(' correct') for line in resp_lines],
                     'Miss': [line.count('incorrect') for line in resp_lines],
                     'RT': [line[line.find('RT')+5:line.rfind(';')].strip() for line in resp_lines],
                     'Timestamp': [float(line[:line.find('.')+4]) for line in resp_lines]})
data['Score'] = [100 if data['Hit'][ix]==1 else 0 for ix in range(len(data))]
data['Condition'] = [line[line.find('condition')+12:line.find('condition')+15] for line in resp_lines]
    
# Calculate ITIs
stim_lines = [line for line in log if line.find('onset =')!=-1]
data['ITI'] = [float(line[line.find('ITI =')+6:line.find('ITI =')+9]) for line in stim_lines]
data.loc[data['Trial']==0,'ITI'] = 0
# data.loc[data['Block']==-1,'ITI'] = 0
            
# Fix Reversals, Block, and missed RTs
for ix in range(len(data)):
#     if data['RT'][ix][0:2]=='-1':#No response
#         data.loc[ix,'RT'] = -1
#         data.loc[ix,'Score'] = 0            # !!! may change depending on version !!!!
#     else:# Real Responses
#         if data['RT'][ix].find(';')!=-1:# shorter number of digits, clip ';'
#             data.loc[ix,'RT'] = float(data['RT'][ix][:data['RT'][ix].find(';')])
#         else:
#             data.loc[ix,'RT'] = float(data['RT'][ix])
    if data['Block'][ix]=='T':# Training
        data.loc[ix,'Block'] = -1
    else:
        data.loc[ix,'Block'] = int(data['Block'][ix])


# In[9]:

# data


# In[10]:

# print(SBJ,' n_trials = ',len(data))
# for cond in data['Condition'].unique():
#     print('==============', cond, '==============')
#     print('Accuracy:', data[data['Condition']==cond]['Hit'].mean())
#     print('--> correct = ', data[data['Condition']==cond]['Hit'].sum(),\
#           '/', len(data[data['Condition']==cond]['Hit']))
#     print('n_reverse:', data[data['Condition']==cond]['Reversal'].sum())


# ### Save Behavioral Data

# In[12]:

behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav_oddball.csv')

data.to_csv(behav_fname,index_label='Total_Trial')


# In[14]:

# pd.set_option('max_rows', 75)
# data[data['Block']==1]


# In[ ]:



