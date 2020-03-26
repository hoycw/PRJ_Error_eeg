function fn_remove_null_channels(SBJ, proc_id)
%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Set up processing and SBJ variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat']);

null_neg = cell(size(SBJ_vars.ch_lab.null));
for null_ix = 1:numel(SBJ_vars.ch_lab.null)
    null_neg{null_ix} = ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.null{null_ix} SBJ_vars.ch_lab.suffix];
end

cfg = [];
cfg.channel = [{'all'}, null_neg];
trials = ft_selectdata(cfg, trials);

clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat'];
save(clean_data_fname, '-v7.3', 'trials', 'cfg_trl', 'ica', 'heog_ics', 'veog_ics', 'eog_trials');
end