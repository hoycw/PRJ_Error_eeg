function EEG02b_ica_rejection(SBJ,proc_id, dorejectvisual)
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

data_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
load(data_fname);
data_cleanname = [SBJ_vars.dirs.preproc SBJ '_clean02a_' proc_id '.mat'];
load(data_cleanname)

%% Import behavioral data
%   Total_Trial,Block,Condition,Hit,RT,Timestamp,Tolerance,Trial,Score,ITI,ITI type
clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav02a_' proc_id '_clean.mat'];
load(clean_bhv_fname);
bad_comp = SBJ_vars.ica_reject;


% IC rejection
cfg = [];
cfg.component = unique([bad_comp, heog_ics, veog_ics]);
cfg.demean = 'no';
clean_trials = ft_rejectcomponent(cfg, ica);

%Repair Bad Channels
cfg = [];
cfg.method = 'average';
cfg.missingchannel = SBJ_vars.ch_lab.bad(:); % not in data (excluded from ica)
cfg.layout = 'biosemi64.lay';
cfgn = [];
cfgn.channel = 'all';
cfgn.layout = 'biosemi64.lay';
cfgn.method = 'template';

cfg.neighbours = ft_prepare_neighbours(cfgn);

clean_trials = ft_channelrepair(cfg, clean_trials);


%% Visual Trial Rejection
if dorejectvisual
    cfg = [];
    cfg.method = 'summary';  % 'summary' for trials+channels; 'channel' for individual trials
    clean_trials = ft_rejectvisual(cfg, clean_trials);
    cfg = [];
    cfg.derivative = 'yes';
    clean_deriv = ft_preprocessing(cfg, clean_trials);
    cfg = [];
    cfg.method = 'summary';
    clean_summ_deriv = ft_rejectvisual(cfg, clean_deriv);
    % Report channels and trials identified above in SBJ_vars, then re-run
else
    bad_final_trials = SBJ_vars.trial_reject_n;
    cfgs = [];
    cfgs.trials = setdiff([1:numel(clean_trials.trial)], bad_final_trials);
    clean_trials = ft_selectdata(cfgs, clean_trials)
end

%% FINAL CHECK
cfg_plot = [];
cfg_plot.viewmode = 'vertical';
out = ft_databrowser(cfg_plot, clean_trials);

bad_final_trials = SBJ_vars.trial_reject_n;
for f_ix = 1:numel(bhv_fields);
     bhv.(bhv_fields{f_ix})(bad_final_trials) = [];
end
%bad_final_trials = input('list trials that still look bad (and add them to bad trials in sbj vars as well!')
%cfgs = [];
%cfgs.trials = setdiff(1:numel(clean_trials.trial),bad_final_trials);
%clean_trials = ft_selectdata(cfgs,clean_trials);
%% Clean up and save data
% Get bad trials and channels from SBJ_vars
%   NOTE: these are indices into the post-raw rejection trial list
%bad_trials = setdiff(SBJ_vars.trial_reject_n+1, exclude_trials')'
%bhv.total_trial is zero indexed so need to add one to get trial number
%for SBJ_vars.trial_reject_n+1
%bad_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad);

% Remove from ICA cleaned data
%cfgs = [];
%cfgs.trials = setdiff(1:numel(clean_trials.trial),bad_trials);
%cfgs.channel = {'all',bad_ch_neg{:}};
%trials = ft_selectdata(cfgs,clean_trials);


    
% Save outputs
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat'];
save(clean_data_fname, '-v7.3', 'clean_trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv');