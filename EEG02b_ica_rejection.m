function EEG02b_ica_rejection(SBJ,proc_id, dorejectvisual)
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
    clean_sum = ft_rejectvisual(cfg, clean_trials);
    cfg = [];
    cfg.derivative = 'yes';
    clean_deriv = ft_preprocessing(cfg, clean_trials);
    cfg = [];
    cfg.method = 'summary';
    clean_summ_deriv = ft_rejectvisual(cfg, clean_deriv);
    cfg_plot.viewmode = 'vertical';
    cfg_plot.ylim = [-15 15];
    ft_databrowser_allowoverlap(cfg_plot, clean_trials);
    % Report channels and trials identified above in SBJ_vars, then re-run
    % these aren't saved to clean trials because that will mess up the
    % indices for the data rejection
else
    fprintf("\nGo run EEG02c please!\n");
end

%% Save outputs
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_clean02b_' proc_id '.mat'];
save(clean_data_fname, '-v7.3', 'clean_trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav02b_' proc_id '_clean.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv', 'bhv_fields');

end