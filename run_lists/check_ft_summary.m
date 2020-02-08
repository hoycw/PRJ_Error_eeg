SBJ = 'EEG23';
block = 1;
proc_id = 'eeg_full_ft';

%%
root_dir = '/Volumes/hoycw_clust/';
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
load([root_dir 'PRJ_Error_eeg/scripts/utils/cfg_plot_eeg.mat']);
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_final.mat']);
% load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02b.mat']);
cfg = [];
cfg.derivative = 'yes';
clean_deriv = ft_preprocessing(cfg, clean_trials);

%%
cfg = [];
cfg.method = 'summary';  % 'summary' for trials+channels; 'channel' for individual trials
clean_sum = ft_rejectvisual(cfg, clean_trials);

%%
cfg = [];
cfg.method = 'summary';
clean_summ_deriv = ft_rejectvisual(cfg, clean_deriv);
