function ODD02c_trial_rejection(SBJ,proc_id, visual)
% This function generates figures for both the ERP stacks and the ICA Plots for the oddball trials.  Also cut out the bad trials (training, RT).
%SBJ = 'EEG#'
%Proc_id = 'egg_full_ft'
%visual = 0 or 1, whether or not to do a final databrowser check of the data

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02b.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_02a.mat']);

%% Eliminate trials
    cfgs = [];
    cfgs.trials = setdiff([1:numel(clean_trials.trial)], SBJ_vars.trial_reject_ix_oddball); % note: trial_reject_ix is the index of the values in clean trials that show up in data_browser.
    %WARNING!! If you do databrowser after doing the steps above and
    %deciding to reject trials in the gooey of lines 44 - 52, your indices
    %will be incorrect.  they need to be the indices of clean_trials in line
    %40. 
    clean_trials = ft_selectdata(cfgs, clean_trials);
%% FINAL CHECK
% Load cfg with plotting parameters
if visual
    cfg_plot.ylim = [-15 15];  %this is different from initial but its in order to see anything unusual
    cfg_plot.viewmode = 'vertical';
    ft_databrowser(cfg_plot, clean_trials);
end

for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(SBJ_vars.trial_reject_ix_oddball) = [];
end

%% Save outputs
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '_final.mat'];
save(clean_data_fname, '-v7.3', 'clean_trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv');


