
# coding: utf-8

# ### This version saves paradigm parameters and does not include some summary variables that can be computed from the main info:
#     reversal, post_err, score_total, hit_total, block_acc

# In[1]:

#get_python().magic(u'matplotlib inline')
import sys 
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as io
import pickle


# In[102]:

SBJ = sys.argv[1]#raw_input('Enter SBJ ID to process:')#'IR63'


# In[103]:
#prj_dir = '/Volumes/hoycw_clust/PRJ_Error_eeg/'
prj_dir = '/Users/sheilasteiner/Desktop/Knight_Lab/PRJ_Error_eeg/'
results_dir = prj_dir+'results/'
fig_type = '.png'
data_dir = prj_dir+'data/'
sbj_dir  = data_dir+SBJ+'/'
# paths = {'Rana': '/Users/colinhoy/Code/PRJ_Error/data/logs/',
#          'Adi': '/Users/colinhoy/Code/PRJ_Error/data/logs/',
#          'IR57': '/Users/colinhoy/Code/PRJ_Error/data/logs/'}
logs = {
        'EP02': 'EEG_pilot2_response_log_20170425110628.txt',
        'EP04': 'Pilot4_response_log_20180322141434.txt',
        'EP05': 'TT_Cyclone_pilot05_response_log_20180426115107.txt',
        'EP06': 'TT_Cyclone_pilot06_response_log_20180426140027.txt',
        'EP07': 'Pilot07_response_log_20181031155619.txt',
        'EP08': 'Pilot08_response_log_20181101084314.txt',
        'EP09': 'Pilot09_response_log_20181127081927.txt',
        'EP10': 'Pilot10_response_log_20181127141339.txt',
        'EP11': 'Pilot11_2_response_log_20181128144407_rm1st3trl.txt',
        'EP14': 'Pilot14_response_log_20190425165756.txt',
	'EP15': 'Pilot15_response_log_20190429183425.txt',
	'EP16': 'Pilot16_response_log_20190430172240.txt',
	'EP17': 'Pilot17_response_log_20190502172924.txt',
	'EP18': 'Pilot18_response_log_20190507110645.txt',
	'EP19': 'Pilot19again_response_log_20190701180728.txt',
        'EEG01': 'eeg01_response_log_20190712151410.txt',
        'EEG02': 'eeg02_response_log_20190715153326.txt',
        'EEG03': 'EEG03_response_log_20190722152821.txt',
        'EEG04': 'eeg04TT_response_log_20190723111922.txt',
        'EEG05': 'eeg05_response_log_20190724153120.txt',
        'EEG06': 'eeg06_response_log_20190730111030.txt',
        'EEG07': 'eeg07startover_response_log_20190802155314.txt',
        'EEG08': 'eeg08_response_log_20190809152412.txt',
        'EEG09': 'eeg09_response_log_20190812175722.txt',
        'EEG10': 'eeg10_response_log_20190815161826.txt',
        'EEG11': 'eeg11_response_log_20190816133923.txt',
        'EEG12': 'eeg12_response_log_20190819143555.txt'
}
#         'EP05': '.txt',# both log files are empty!
# logs = {'Rana_1.6': 'Rana2_response_log_20170321103129_DATA.txt',
#         'Adi_1.7': 'adi_response_log_20170321153641.txt',
#         'CP22': '222_response_log_20170609140407.txt',
#         'CP23': '223_response_log_20170930123015.txt',
#         'CP24': '224_response_log_20171206121023.txt',
#         'CP242': 'cp24_2_response_log_20171209120902.txt',
#         'CP25': 'CP25_response_log_20180811082953.txt',
#         'IR57': '857_response_log_20170322112243_CWHedit.txt',#CWHedit added n_training, n_examples lines
#         'IR60': '60_response_log_20170613100307.txt',
#         'IR62': 'ir62_response_log_20170713124719.txt',
#         'IR63': 'IR63_response_log_20170921095757.txt',
#         'IR65': '865_response_log_20171207130759.txt',
#         'IR66': 'ir66_response_log_20171219124409.txt',
#         'IR67': 'ir673_response_log_20180124103600.txt',
#         'IR68': 'IR68_response_log_20180124140950.txt',
#         'IR69': '869p2_response_log_20180211111609.txt',#this is 2nd run, another one before!
#         'IR71': '71_response_log_20180221115510.txt',
#         'IR72': 'Ir72_response_log_20180325133520.txt',#2 log files, this has easy, other has hard blocks
#         'IR74': 'ir742_response_log_20180327170333.txt',
#         'IR75': 'IR75_response_log_20180531221813.txt',
#         'IR76': 'IR76_response_log_20180603181356.txt',
#         'IR77': '877_response_log_20180620121202.txt',#1/3 logs, has 2 easy, 8773 has 2 hards
#         'IR78': 'IR78_response_log_20180628052715.txt',
#         'IR79': 'IR79_response_log_20180710112314.txt',
#         'IR82': 'IR82_response_log_20180928162311.txt',
#         'IR84': 'IR84_response_log_20181025094454.txt',
#         'BP1': 'Pilot1_2_response_log_20170412131644.txt',
#         'BP2': 'pilot2_response_log_20170412140615.txt',
#         'BP3': 'pilot3_response_log_20170413110228.txt',
#         'BP4': 'Pilot4_2_response_log_20170418140941.txt',
#         'BP5': 'colin_real_response_log_20170412103113.txt',
#         'BP6': 'pilot_adi_response_log_20170414122257.txt',
#         'BP7': 'pilot_Rana_response_log_20170415155844.txt',
#         'BP8': 'Giao_response_log_20170419161340.txt',
#         'BP9': 'Sundberg_response_log_20170419150222.txt',
#         'colin_vec': 'colin_circle_wVec_response_log_20171222141248.txt',
#         'colin_novec': 'colin_noVec_response_log_20171222142110.txt'
#        }


# ### Load SBJ Data

# In[104]:

log_filename = os.path.join(sbj_dir,'00_raw',logs[SBJ])

log_file = open(log_filename,'r')
log = log_file.readlines()
log_file.close()


# ### Process Version-Specific Parameters

# In[105]:

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
        prdm['ITIs'] = [float(string) for string in line[line.find('[')+1:line.find(']')].split(',')]
        ITI_bounds = np.mean([prdm['ITIs'][:-1], prdm['ITIs'][1:]],0)
        
    # Tolerance limits/clamps
    if line.find('tolerance_lim')!=-1:
        prdm['tol_lim'] = [float(string) for string in line[line.find('[')+1:line.find(']')].split(',')]
        
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

prdm['trl_len'] = prdm['target'] + prdm['fb_delay'] + prdm['fb']


# ### Save paradigm parameters

# In[106]:

# Python readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars.pkl')
with open(prdm_fname, 'wb') as f:
    pickle.dump(prdm, f, pickle.HIGHEST_PROTOCOL)
# MATLAB Readable
prdm_fname = os.path.join(sbj_dir,'03_events',SBJ+'_prdm_vars.mat')
io.savemat(prdm_fname,prdm)


# In[107]:

print ('paradigm: ', prdm['prdm_name'], ' v', prdm['prdm_version'])
print
print ('interval: ', prdm['target'], 's')
print ('feedback_delay: ', prdm['fb_delay'], 's')
print ('feedback duration: ', prdm['fb'], 's')
print ('total trial length: ', prdm['trl_len'], 's')
print
print ('n_blocks: ', prdm['n_blocks'])
print ('n_trials/block: ', prdm['n_trials'])
print ('n_full_vis_examples: ', prdm['n_examples'])
print ('n_training/condition: ', prdm['n_training'])
print
# ITI_bounds = [np.mean(a,b) for a, b in zip(ITIs[:-1],ITIs[1:])]
print ('ITIs:',prdm['ITIs'], ITI_bounds)
print ('tolerance_lim:', prdm['tol_lim'])


# ### Extract Trial Info

# In[108]:

resp_lines = [line for line in log if line.find('Outcome=')!=-1]
data = pd.DataFrame({'Block': [line[line.find('B')+1] for line in resp_lines],
                     'Trial': [int(line[line.find('_T')+2:line.find(':')]) for line in resp_lines],
                     'Condition': [line[line.find('_type')+8:line.find('_type')+12] for line in resp_lines],
                     'Hit': [line.count('WIN') for line in resp_lines],
                     'RT': [line[line.find('RT')+5:line.find('RT')+5+13].strip() for line in resp_lines],
                     'Tolerance': [float(line[line.find('tol')+12:line.find('\n')]) for line in resp_lines],
                     'Timestamp': [float(line[:line.find('.')+4]) for line in resp_lines]})
data['Score'] = [100 if data['Hit'][ix]==1 else -100 for ix in range(len(data))]
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
            
# Fix Reversals, Block, and missed RTs
for ix in range(len(data)):
    if data['RT'][ix][0:2]=='-1':#No response
        data.loc[ix,'RT'] = -1
        data.loc[ix,'Score'] = 0            # !!! may change depending on version !!!!
    else:# Real Responses
        if data['RT'][ix].find(';')!=-1:# shorter number of digits, clip ';'
            data.loc[ix,'RT'] = float(data['RT'][ix][:data['RT'][ix].find(';')])
        else:
            data.loc[ix,'RT'] = float(data['RT'][ix])
    if data['Block'][ix]=='T':# Training
        data.loc[ix,'Block'] = -1
    else:
        data.loc[ix,'Block'] = int(data['Block'][ix])


# In[109]:

# print(SBJ,' n_trials = ',len(data))
# for cond in data['Condition'].unique():
#     print('==============', cond, '==============')
#     print('Accuracy:', data[data['Condition']==cond]['Hit'].mean())
#     print('--> correct = ', data[data['Condition']==cond]['Hit'].sum(),\
#           '/', len(data[data['Condition']==cond]['Hit']))
#     print('n_reverse:', data[data['Condition']==cond]['Reversal'].sum())


# ### Save Behavioral Data

# In[110]:

behav_fname = os.path.join(sbj_dir,'03_events',SBJ+'_behav.csv')

data.to_csv(behav_fname,index_label='Total_Trial')


# In[111]:

# pd.set_option('max_rows', 75)
# data[data['Block']==1]


# In[ ]:



