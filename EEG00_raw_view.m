function EEG00_raw_view(SBJ)
%% Only for viewing raw data and marking epochs to toss

%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath(genpath([root_dir 'PRJ_Error_eeg/scripts/']));
addpath(ft_dir);
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

cfg=[];
cfg.dataset  = SBJ_vars.dirs.raw_filename;
cfg.lpfilter = 'no';
cfg.hpfilter = 'no';
cfg.bpfilter = 'yes';
cfg.bpfreq   = [0.1 200];
raw = ft_preprocessing(cfg);

%% Plot and mark bad epochs
load([root_dir 'PRJ_Error_eeg/scripts/utils/cfg_plot_eeg.mat']);
browsed_raw = ft_databrowser(cfg_plot, raw);
bad_epochs = browsed_raw.artfctdef.visual.artifact;

%% Save out bad epochs
out_fname = [SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat'];
save(out_fname, 'bad_epochs');

end