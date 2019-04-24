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
clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav02a_' proc_id '_clean.mat'];
load(clean_bhv_fname);

%% IC rejection
cfg = [];
cfg.component = unique([SBJ_vars.ica_reject, heog_ics, veog_ics]);
cfg.demean = 'no';
clean_trials = ft_rejectcomponent(cfg, ica);

%% Repair Bad Channels
cfg = [];
cfg.method         = 'average';
cfg.missingchannel = SBJ_vars.ch_lab.bad(:); % not in data (excluded from ica)
cfg.layout         = 'biosemi64.lay';

cfgn = [];
cfgn.channel = 'all';
cfgn.layout  = 'biosemi64.lay';
cfgn.method  = 'template';
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
    cfgs = [];
    cfgs.trials = setdiff([1:numel(clean_trials.trial)], SBJ_vars.trial_reject_n);
    clean_trials = ft_selectdata(cfgs, clean_trials);
end

%% FINAL CHECK
cfg_plot = [];
cfg_plot.viewmode = 'vertical';
out = ft_databrowser(cfg_plot, clean_trials);   %SHEILA: do you use the out here? what are you checking here?

% SHEILA: you're confusing trial_n and ix again here- since some trials are
% already gone from the bhv struct, using _n here will be wrong...
for f_ix = 1:numel(bhv_fields);
     bhv.(bhv_fields{f_ix})(SBJ_vars.trial_reject_n) = [];
end

%% Save outputs
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat'];
save(clean_data_fname, '-v7.3', 'clean_trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv');

end