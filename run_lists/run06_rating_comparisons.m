%% Behavior and Pre-Processing for Sequential PE Initial Submission
% Written 8/14/2020 by Colin Hoy
% 	Fig. 1A: BHV02_plot_RT_hist
%   Fig. 1B: BHV00_01_prelim_analysis.py
%   Fig. 1C: BHV03_group_accuracy_plots_TT

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
SBJ_id = 'test_ratings';
SBJs = fn_load_SBJ_list(SBJ_id);

%% Pre-Processing
% Rating behavior data pre-processing
%   RATE01_behavior_rejection.m
%       -Cut out training and bad RT trials from behavior struct
%       -save final preprocessed data and behavior

proc_id                  = 'eeg_full_ft';
fig_vis                  = 'on';
% odd_proc_id = 'odd_full_ft';  % for Oddball task preprocessing

% NOTE: This was never run cleanly, but SBJ by SBJ and occasionally in a loop
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Toss training and bad RTs
    RATE01_behavior_rejection(SBJs{s}, proc_id);
    
    clear SBJ_vars
end

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
%***BHV03_plot_group_accuracy_TT.m: group accuracy bar plot for Fig. 1C

% Fig. 1B: Call BHV00 and BHV01
%   usually run on command line, potentially commenting out BHV00/01 in .sh
bhv_cmd = ['bash run_BHV00_01_prelim_analysis.sh << ' SBJ_id];
% system(bhv_cmd);

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_ftype = 'png';

% Fig. 1A: Example RT Histogram
SBJ       = 'colin_test1';
BHV02_plot_RT_hist(SBJ,proc_id,save_fig,'fig_ftype',fig_ftype);

% Fig. 1C: Group Accuracy
conditions = 'Dif';
% BHV03_plot_group_accuracy_TT(SBJ_id, conditions, fig_ftype);

%% Survey Data Analysis
%   Test neutral trial survey responses for valence
plot_fig  = 1;
save_fig  = 0;
fig_ftype = 'png';

% Run only for 24 good 'EEG' SBJs, which filled out final version of survey
% BHV04_survey_grp_stats('goodEEG', plot_fig, save_fig, fig_ftype);

%% Single SBJ RL Model
proc_id   = 'eeg_full_ft';
stat_ids  = {'ERPEs_DifFB_lme_st05t5'};%'uRPE_Neg_lme_st05t5'};%'ML_Neg_lme_st05t5'};%'ERPEsL_all_lme_st05t5'};
% Alternative (worse) models: 'RSVPE_all_lme_mn1FRN','SML_all_lme_mn1FRN','VML_all_lme_mn1FRN'
fig_vis   = 'on';
save_fig  = 1;
fig_ftype = 'png';

for s = 1:numel(SBJs)
    for st_ix = 1:numel(stat_ids)
        % Run model
%         SBJ04a_RL_model_ratings(SBJs{s},proc_id,stat_ids{st_ix});
        
        % Fig. 1D: Plot model fit to tolerance and outcomes/accuracy
        SBJ04b_BHV_RL_model_rating_plot(SBJs{s},proc_id,stat_ids{st_ix},...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    % close all;
end

% Fig. 1D inset: Plot group model fits (overlapping sigmoids without tolerance)
% SBJ04b_BHV_RL_model_plot_grp(SBJ_id,proc_id,stat_id,...
%     'fig_vis',fig_vis,'fig_ftype',fig_ftype);

% Sup. Fig. 1A: Plot model predicitons by condition across group
plt_id    = 'line_cond';
for st_ix = 1:numel(stat_ids)
%     SBJ04b_plot_model_predictions(SBJ_id,proc_id,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

