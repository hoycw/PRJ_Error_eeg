
# coding: utf-8

# ### This version saves paradigm parameters and does not include some summary variables that can be computed from the main info:
#     reversal, post_err, score_total, hit_total, block_acc

# In[1]:

#%matplotlib inline
import sys 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as io
import pickle


# In[2]:

SBJ = sys.argv[1]#raw_input('Enter SBJ ID to process:')#'EEG01'


# In[3]:
prj_dir = '/Users/sheilasteiner/Desktop/Knight_Lab/PRJ_Error_eeg/'
#prj_dir = '/Volumes/hoycw_clust/PRJ_Error_eeg/'
results_dir = prj_dir+'results/'
fig_type = '.png'
data_dir = prj_dir+'data/'
sbj_dir  = data_dir+SBJ+'/'

logs = {}
with open(data_dir+'TT_behav_log_list.txt') as f:
    for line in f:
        (sbj, log_name) = line.split(',')
        logs[sbj] = log_name.replace('\n','')


# ### Load SBJ Data

# In[4]:

log_fname = os.path.join(sbj_dir,'00_raw',logs[SBJ])

log_file = open(log_fname,'r')
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

    # Timing variables
    if line.find('interval_dur =')!=-1:
        prdm['target'] = float(line[line.find('= ')+2:line.find('\n')])
    if line.find('feedback_delay =')!=-1:
        prdm['fb_delay'] = float(line[line.find('= ')+2:line.find('\n')])
    if line.find('feedback_dur =')!=-1:
        prdm['fb'] = float(line[line.find('= ')+2:line.find('\n')])
    
    # ITIs and boundaries between them
    if line.find('ITIs')!=-1:
        prdm['ITIs'] = [float(string)                     for string in line[line.find('[')+1:line.find(']')].split(',')]
        ITI_bounds = np.mean([prdm['ITIs'][:-1], prdm['ITIs'][1:]],0)
        
    # Tolerance limits/clamps
    if line.find('tolerance_lim')!=-1:
        prdm['tol_lim'] = [float(string)                          for string in line[line.find('[')+1:line.find(']')].split(',')]
        
    # Trial count variables
    if line.find('n_blocks')!=-1:
        prdm['n_blocks'] = int(line[line.find('=')+2:])
    if line.find('n_trials')!=-1:
        prdm['n_trials'] = int(line[line.find('=')+2:])
    if line.find('n_examples')!=-1:
        prdm['n_examples'] = int(line[line.find('=')+2:])
    elif line.find('n_fullvis')!=-1:
        prdm['n_examples'] = int(line[line.find('=')+2:])
    if line.find('n_training')!=-1:
        prdm['n_training'] = int(line[line.find('=')+2:])

# Add missing items from early log files
if 'prdm_name' not in prdm:
    prdm['prdm_name'] = 'not_logged'
if 'prdm_version' not in prdm:
    prdm['prdm_version'] = '<1.8.5'
if 'n_examples' not in prdm:
    prdm['n_examples'] = int(-1)
if 'n_training' not in prdm:
    prdm['n_training'] = int(-1)

prdm['trl_len'] = prdm['target']+            prdm['fb_delay']+prdm['fb']


# ### Save paradigm parameters

# In[6]:

# Python readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars.pkl')
with open(prdm_fname, 'wb') as f:
    pickle.dump(prdm, f, pickle.HIGHEST_PROTOCOL)
# MATLAB Readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars.mat')
io.savemat(prdm_fname,prdm)


# In[7]:

print 'paradigm: ', prdm['prdm_name'], ' v', prdm['prdm_version']
print
print 'interval: ', prdm['target'], 's'
print 'feedback_delay: ', prdm['fb_delay'], 's'
print 'feedback duration: ', prdm['fb'], 's'
print 'total trial length: ', prdm['trl_len'], 's'
print
print 'n_blocks: ', prdm['n_blocks']
print 'n_trials/block: ', prdm['n_trials']
print 'n_full_vis_examples: ', prdm['n_examples']
print 'n_training/condition: ', prdm['n_training']
print
# ITI_bounds = [np.mean(a,b) for a, b in zip(ITIs[:-1],ITIs[1:])]
print 'ITIs:',prdm['ITIs'], ITI_bounds
print 'tolerance_lim:', prdm['tol_lim']


# ### Extract Trial Info

# In[8]:

resp_lines = [line for line in log if line.find('Outcome=')!=-1]
data = pd.DataFrame({'Block': [line[line.find('B')+1] for line in resp_lines],
                     'Trial': [int(line[line.find('_T')+2:line.find(':')]) for line in resp_lines],
                     'Feedback': ['W' if line.count('WIN')>0 else \
                                  'L' if line.count('LOSE')>0 else \
                                  'S' for line in resp_lines],
                     'RT': [line[line.find('RT')+5:line.find('RT')+5+13].strip() for line in resp_lines],
                     'Tolerance': [float(line[line.find('tol')+12:line.find('\n')]) for line in resp_lines],
                     'Timestamp': [float(line[:line.find('.')+4]) for line in resp_lines]})

# Fix Reversals, Block, and missed RTs
for ix in range(len(data)):
    # Fix RTs
    if data['RT'][ix][0:2]=='-1':#No response
        data.loc[ix,'RT'] = -1
        #data.loc[ix,'Score'] = 0            # !!! may change depending on version !!!!
    else:# Real Responses
        if data['RT'][ix].find(';')!=-1:# shorter number of digits, clip ';'
            data.loc[ix,'RT'] = float(data['RT'][ix][:data['RT'][ix].find(';')])
        else:
            data.loc[ix,'RT'] = float(data['RT'][ix])
    
    # Fix Block coding
    if data['Block'][ix]=='T':# Training
        data.loc[ix,'Block'] = -1
    else:
        data.loc[ix,'Block'] = int(data['Block'][ix])

data['Hit'] = [1 if abs(data['RT'][ix]-prdm['target']) <= data['Tolerance'][ix] else 0 for ix in range(len(data))]
data['Score'] = [0 if data['Feedback'][ix]=='S' or data['RT'][ix]== -1 else                  100 if data['Feedback'][ix]=='W' else -100 for ix in range(len(data))]

# Mark trials with bad feedback (apparently my logic has error up to 10 ms... damn it!)
data['bad_fb'] = [False if data['Feedback'][ix]=='S' or data['RT'][ix]==-1 or                   (data['Feedback'][ix]=='W' and data['Hit'][ix]==1) or                   (data['Feedback'][ix]=='L' and data['Hit'][ix]==0)                   else True for ix in range(len(data))]

# Add condition based on logging version
if prdm['prdm_version'][0]=='2':
    data['Condition'] = [line[line.find('condition')+12:line.find('condition')+16] for line in resp_lines]
else:
    data['Condition'] = [line[line.find('_type')+8:line.find('_type')+12] for line in resp_lines]
    
# Calculate ITIs
data['ITI'] = [data['Timestamp'][ix]-data['Timestamp'][ix-1]-prdm['trl_len'] if ix!=0 else 0                for ix in range(len(data))]
data.loc[data['Trial']==0,'ITI'] = 0
data.loc[data['Block']==-1,'ITI'] = 0
# Match real ITIs to ITI categories
ITI_bin_edges = np.insert(ITI_bounds, ITI_bounds.shape[0], prdm['ITIs'][-1]+0.1)
ITI_bin_edges = np.insert(ITI_bin_edges, 0, 0)
data['ITI type'] = [prdm['ITIs'][np.argmax(np.histogram(data['ITI'][ix],bins=ITI_bin_edges)[0])]        if data['ITI'][ix]!=0 else 0 for ix in range(len(data))]
# if len(prdm['ITIs'])==4:    # target_time v1.8.5+
# elif len(prdm['ITIs'])==3:  # target_time v1.8.4 and below
# else:               # Errors for anything besides len(ITIs)==3,4
#     assert len(prdm['ITIs'])==4
            
if any(data['bad_fb']):
    bad_ix = [ix for ix in range(len(data)) if data['bad_fb'][ix]]
    tmp = data.ix[bad_ix]
    tmp['error'] = abs(tmp['RT']-prdm['target'])
    tmp['error-tol'] = tmp['Tolerance']-tmp['error']
    print 'WARNING!!! Bad feedback found on {0} trials!'.format(len(bad_ix))
    print 'Max error-tolerance = {0}, mean = {1}'.format(np.max(tmp['error-tol']),np.mean(np.abs(tmp['error-tol'])))
#   tmp.ix[:,{'Tolerance','error','WIN','Hit','err-tol'}]


# In[9]:

data


# In[10]:

# print(SBJ,' n_trials = ',len(data))
# for cond in data['Condition'].unique():
#     print('==============', cond, '==============')
#     print('Accuracy:', data[data['Condition']==cond]['Hit'].mean())
#     print('--> correct = ', data[data['Condition']==cond]['Hit'].sum(),\
#           '/', len(data[data['Condition']==cond]['Hit']))
#     print('n_reverse:', data[data['Condition']==cond]['Reversal'].sum())


# ### Save Behavioral Data

# In[11]:

behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav.csv')

data.to_csv(behav_fname,index_label='Total_Trial')


# In[131]:

# pd.set_option('max_rows', 75)
# data[data['Block']==1]


# In[ ]:



