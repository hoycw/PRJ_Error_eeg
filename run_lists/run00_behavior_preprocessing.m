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
% EEG data pre-processing
%   SBJ00_raw_view.m: view raw data to mark worst epochs
%       -load, filter, downsample, plot PSDs
%       -plot data to mark and save worst epochs to exclude from ICA
%   SBJ01_preproc.m: preprocess the data
%       -import, preprocess, and concatenate data blocks
%       -NaN out bad epochs from SBJ00
%       -run ICA on data with NaNs to avoid artifacts dominating ICs
%       -save data without NaNs (to allow filtering on whole time series) and ICA results
%       NOTE: see comments in scripts/run_lists/rerun_preprocessing/resave_SBJ01_data_without_NaNs.m
%   SBJ02a_artifact_rejection.m
%       -Segment trials
%       -Cut out raw bad_epoch, training, and bad RT trials from data and behavior struct
%       -Identify EOG ICs via correlation
%       -Generate figures for quality checks on ERPs and ICA components.
%   SBJ02b_ica_rejection.m
%   SBJ02c_trail_rejection.m

proc_id                  = 'eeg_full_ft';
view_previous_bad_epochs = 1;
generate_figs            = 1;
fig_vis                  = 'on';
clear_prev_QA_plots      = 0;
% odd_proc_id = 'odd_full_ft';  % for Oddball task preprocessing

% NOTE: This was never run cleanly, but SBJ by SBJ and occasionally in a loop
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % View raw data
    for b_ix = 1:numel(SBJ_vars.block_name)
        SBJ00_raw_view(SBJs{s},view_previous_bad_epochs,proc_id,b_ix);
    end
    % Preprocess data and run ICA
    SBJ01_preproc(SBJs{s},proc_id);
    % 
    SBJ02a_artifact_rejection(SBJs{s}, proc_id, generate_figs, fig_vis, clear_prev_QA_plots);
    % 
    SBJ02b_ica_rejection(SBJs{s}, proc_id, reject_visual);
    % 
    SBJ02c_trial_rejection(SBJs{s}, proc_id, plot_final_check);
    
    clear SBJ_vars b_ix
end