%% Behavior and Pre-Processing for Sequential PE Initial Submission
% Written 8/14/2020 by Colin Hoy

%% Set Up
% Running locally on Colin's MacBookPro; OS 10.13.6; MATLAB R2017b
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Subject Lists
% root_dir/PRJ_Error_EEG/scripts/SBJ_lists contains text files with SBJ IDs
%   -good1.sbj: cohort 1
%   -good2.sbj: cohort 2
%   -goodall.sbj: cohort 1+2
%   -goodEEG1.sbj: cohort 1 (only with oddball)
%   -goodEEG2.sbj: cohort 2 (only with oddball; same as good2.sbj)
%   -goodEEG.sbj: cohort 1+2 (only with oddball)
SBJ_id = 'goodall';

%% Behavior (Fig. 1B, 1C)
% Behavioral pre-processing and analysis
%   run_BHV00_01_prelim_analysis.sh: pulls SBJ_list to automatically run BHV00 and BHV01 for those SBJs
%   runs in PRJ_Error_py2.7 conda env with python 2.7
%       -colin_PRJ_error_py2.7.yml: contains conda env package list
%   BHV00_extract.ipynb (and corresponding .py): 
%       -root_dir/PRJ_Error_eeg/data/TT_bheav_log_list.txt: text file containing names of behavioral logs per SBJ
%       -read, save, print paradigm parameters (python and MATLAB readable)
%       -read, pre-process, QA check, and save behavioral data as CSV
%   BHV01_prelim_analysis.ipynb (and corresponding .py): 
%       -load behavioral data
%       -correct idiosyncracies (block numbers, extra training, missing trials, etc.)
%       -compute accuracy
%***    -plot Fig. 1B SBJ performance over time  
%       -plot QA checks: ITI histogram
%       -plot behavioral checks: RT histogram, RT histogram by ITI, RT
%       after long vs. short RTs, accuracy by ITI and condition
%***BHV02_plot_RT_hist.m: plot RT histogram for Fig. 1A
%***BHV03_plot_group_accuracy_TT.m

% Fig. 1B: Call BHV00 and BHV01
%   usually run on command line, potentially commenting out BHV00/01 in .sh
% bhv_cmd = ['bash run_BHV00_01_prelim_analysis.sh << ' SBJ_id];
% system(bhv_cmd);

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_ftype = 'svg';

% Fig. 1A: Example RT Histogram
SBJ       = 'EEG13';
BHV02_plot_RT_hist(SBJ,proc_id,save_fig,'fig_ftype',fig_ftype);

% Fig. 1C: Group Accuracy
conditions = 'Dif';
BHV03_group_accuracy_plots_TT(SBJ_id, conditions, fig_ftype);

%% Pre-Processing
