function SBJ02c_trial_rejection(SBJ, proc_id, plot_final_check)
% This function generates figures for both the ERP stacks and the ICA Plots.  Also cut out the bad trials (training, RT).
% INPUTS:
%   SBJ [str] - name of the SBJ
%   proc_id [str] - name of the preprocessing pipeline parameters (e.g., 'egg_full_ft')
%   plot_final_check [0/1] - binary flag to run final databrowser check of the data
% OUTPUTS:
%   clean_trials [FT struct] - final preprocessed data after rejecting trials
%   bhv [struct] - final preprocessed behavior after rejecting trials

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
% NOTE: SBJ_vars.trial_reject_ix is the index of the values in clean_trials
%   as shown in ft_databrowser, i.e., after rejecting bad RTs, training, and bad_epochs
% WARNING!! Multiple iterations rejecting trials using ft_databrowser "gooey"
%   will result in incorrect indices. All indices need to be into
%   clean_trials from SBJ02b_ica_rejection.
cfgs = [];
cfgs.trials = setdiff([1:numel(clean_trials.trial)], SBJ_vars.trial_reject_ix);
clean_trials = ft_selectdata(cfgs, clean_trials);

%% FINAL CHECK
% Load cfg with plotting parameters
if plot_final_check
    cfg_plot.ylim = [-15 15];  %this is different from initial but its in order to see anything unusual
    cfg_plot.viewmode = 'vertical';
    ft_databrowser(cfg_plot, clean_trials);
end

% Remove bad trials from behavior struct
for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(SBJ_vars.trial_reject_ix) = [];
end

%% Save outputs
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '_final.mat'];
save(clean_data_fname, '-v7.3', 'clean_trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv');

end
