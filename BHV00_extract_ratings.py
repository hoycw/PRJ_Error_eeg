
# coding: utf-8

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


# In[28]:

SBJ = sys.argv[1]#raw_input('Enter SBJ ID to process:')#'EEG01'


# In[29]:

prj_dir = '/Volumes/hoycw_clust/PRJ_Error_eeg/'
results_dir = prj_dir+'results/'
fig_type = '.png'
data_dir = prj_dir+'data/'
sbj_dir  = data_dir+SBJ+'/'

# Get log names
logs = {}
with open(data_dir+'TT_ratings_behav_log_list.txt') as f:
    for line in f:
        (sbj, log_name) = line.split(',')
        logs[sbj] = log_name.replace('\n','')


# ### Load SBJ Data

# In[30]:

log_fname = os.path.join(sbj_dir,'00_raw',logs[SBJ])

log_file = open(log_fname,'r')
log = log_file.readlines()
log_file.close()


# ### Process Version-Specific Parameters

# In[31]:

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
    if line.find('bad_fb_tolerance =')!=-1:
        prdm['bad_fb_tolerance'] = float(line[line.find('= ')+2:line.find('\n')])
    if line.find('rating_delay =')!=-1:
        prdm['rating_delay'] = float(line[line.find('= ')+2:line.find('\n')])
    if line.find('max_rating_time =')!=-1:
        prdm['max_rating_time'] = float(line[line.find('= ')+2:line.find('\n')])
    
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
prdm['rating_trl_len'] = prdm['target']+prdm['rating_delay']+            prdm['fb_delay']+prdm['fb']


# ### Save paradigm parameters

# In[32]:

# Python readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars.pkl')
with open(prdm_fname, 'wb') as f:
    pickle.dump(prdm, f, pickle.HIGHEST_PROTOCOL)
# MATLAB Readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars.mat')
io.savemat(prdm_fname,prdm)


# In[33]:

## print 'paradigm: ', prdm['prdm_name'], ' v', prdm['prdm_version']
## print
## print 'interval: ', prdm['target'], 's'
## print 'feedback_delay: ', prdm['fb_delay'], 's'
## print 'rating_delay: ', prdm['rating_delay'], 's'
## print 'max_rating_time: ', prdm['max_rating_time'], 's'
## print 'feedback duration: ', prdm['fb'], 's'
## print 'total trial length: ', prdm['trl_len'], 's'
## print
## print 'n_blocks: ', prdm['n_blocks']
## print 'n_trials/block: ', prdm['n_trials']
## print 'n_full_vis_examples: ', prdm['n_examples']
## print 'n_training/condition: ', prdm['n_training']
## print
## # ITI_bounds = [np.mean(a,b) for a, b in zip(ITIs[:-1],ITIs[1:])]
## print 'ITIs:',prdm['ITIs'], ITI_bounds
## print 'tolerance_lim:', prdm['tol_lim']


# ### Extract Trial Info

# In[48]:

# Extract rating data
if SBJ=='colin_test1':
    rating_lines = [line for line in log if line.find('Rating=')!=-1]
    ratings = [int(line[line.find('Rating=')+7:line.find('Rating=')+9].replace(';','')) for line in rating_lines]
else:
    rating_lines = [line for line in log if line.find('Rating =')!=-1]
    ratings = [int(line[line.find('Rating =')+9:line.find('Rating =')+12].replace(';','')) for line in rating_lines]

rating_rts = [float(line[line.find('rating_RT =')+12:line.find('rating_RT =')+16].replace(';','')) for line in rating_lines]

# get block_n and trial_n assuming no ratings in BT (training)
rating_blocks = [int(line[line.find('B')+1]) for line in rating_lines]
rating_trial_n = [int(line[line.find('_T')+2:line.find(':')]) for line in rating_lines]


# In[23]:

resp_lines = [line for line in log if line.find('Outcome=')!=-1]
data = pd.DataFrame({'Block': [line[line.find('B')+1] for line in resp_lines],
                     'Trial': [int(line[line.find('_T')+2:line.find(':')]) for line in resp_lines],
                     'Feedback': ['W' if line.count('WIN')>0 else \
                                  'L' if line.count('LOSE')>0 else \
                                  'S' for line in resp_lines],
                     'RT': [line[line.find('RT')+5:line.find('RT')+5+13].strip() for line in resp_lines],
                     'Tolerance': [float(line[line.find('tol')+12:line.find('\n')]) for line in resp_lines],
                     'Timestamp': [float(line[:line.find('.')+4]) for line in resp_lines],
                     'Rating': [-1 for line in resp_lines],
                     'Rating_RT': [-1 for line in resp_lines]
                    })

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

# Add ratings and rating_RTs
rating_ix = 0
for ix in range(len(data)):
    if data.loc[ix,'Block']==rating_blocks[rating_ix] and data.loc[ix,'Trial']==rating_trial_n[rating_ix]:
        data.loc[ix,'Rating'] = ratings[rating_ix]
        data.loc[ix,'Rating_RT'] = rating_rts[rating_ix]
        rating_ix += 1
    if rating_ix == len(ratings):
        break


data['Hit'] = [1 if abs(data['RT'][ix]-prdm['target']) <= data['Tolerance'][ix] else 0 for ix in range(len(data))]
data['Score'] = [0 if data['Feedback'][ix]=='S' or data['RT'][ix]== -1 else                  100 if data['Feedback'][ix]=='W' else -100 for ix in range(len(data))]

# Mark trials with bad feedback (apparently my logic has error up to 10 ms... damn it!)
data['bad_fb'] = [False if data['Feedback'][ix]=='S' or data['RT'][ix]==-1 or                   (data['Feedback'][ix]=='W' and data['Hit'][ix]==1) or                   (data['Feedback'][ix]=='L' and data['Hit'][ix]==0)                   else True for ix in range(len(data))]

# Add condition based on logging version
if prdm['prdm_version'][0]=='2' or prdm['prdm_version'] in ['3.0','3.1']:
    data['Condition'] = [line[line.find('condition')+12:line.find('condition')+16] for line in resp_lines]
else:
    data['Condition'] = [line[line.find('_type')+8:line.find('_type')+12] for line in resp_lines]
    
# Calculate ITIs starting on second trial, accounting for ratings
data['ITI'] = [0 if ix==0 else                data['Timestamp'][ix]-data['Timestamp'][ix-1]-prdm['trl_len'] if data['Rating'][ix]==-1                else data['Timestamp'][ix]-data['Timestamp'][ix-1]-prdm['rating_trl_len']-data['Rating_RT'][ix]                for ix in range(len(data))]

# first test run was messed up, 2x rating_delay on all non-rating trials too
if SBJ=='colin_test1':
    for ix in range(len(data)):
        if data['Rating'][ix]==-1:
            data.loc[ix,'ITI'] = data.loc[ix,'ITI']-prdm['rating_delay']

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
            
# Print stats on bad feedback
if any(data['bad_fb']):
    bad_ix = [ix for ix in range(len(data)) if data['bad_fb'][ix]]
    tmp = data.ix[bad_ix]
    tmp['error'] = abs(tmp['RT']-prdm['target'])
    tmp['error-tol'] = tmp['Tolerance']-tmp['error']
    print 'WARNING!!! Bad feedback found on {0} trials!'.format(len(bad_ix))
    print 'Max error-tolerance = {0}, mean = {1}'.format(np.max(tmp['error-tol']),np.mean(np.abs(tmp['error-tol'])))
    tmp.ix[:,{'Tolerance','error','WIN','Hit','err-tol'}]


# In[24]:

# print(SBJ,' n_trials = ',len(data))
# for cond in data['Condition'].unique():
#     print('==============', cond, '==============')
#     print('Accuracy:', data[data['Condition']==cond]['Hit'].mean())
#     print('--> correct = ', data[data['Condition']==cond]['Hit'].sum(),\
#           '/', len(data[data['Condition']==cond]['Hit']))
#     print('n_reverse:', data[data['Condition']==cond]['Reversal'].sum())


# ### Save Behavioral Data

# In[25]:

behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav.csv')

data.to_csv(behav_fname,index_label='Total_Trial')


# In[48]:

# Check ITI computation and timing
# ix = 42
# if data['Rating'][ix]==-1:
#     print 'no rating: ', data['Rating'][ix], data['Rating_RT'][ix], 'last_t:', data['Timestamp'][ix-1], \
#         'curr_t:', data['Timestamp'][ix], prdm['trl_len'], data['Timestamp'][ix]-data['Timestamp'][ix-1]-prdm['trl_len']
# else:
#     print 'ys rating: ', data['Rating'][ix], data['Rating_RT'][ix], 'last_t:', data['Timestamp'][ix-1], 'curr_t:', data['Timestamp'][ix], \
#         data['Timestamp'][ix]-data['Timestamp'][ix-1]-prdm['rating_trl_len']-data['Rating_RT'][ix]


# In[15]:

## pd.set_option('max_rows', 75)
## data[data['Block']==5]


# In[ ]:



