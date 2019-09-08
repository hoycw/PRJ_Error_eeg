function oddball_stats_plot(SBJ, proc_id, conditions)
%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/', ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load preprocessed data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_oddball_' proc_id '_clean.mat'];
load(clean_bhv_fname);
data_cleanname = [SBJ_vars.dirs.preproc SBJ '_clean_oddball_' proc_id '.mat'];
load(data_cleanname)

%find index of conditions of interest
std_cdx = find(clean_oddball.trialinfo == 1);
tar_cdx = find(clean_oddball.trialinfo == 2);
odd_cdx = find(clean_oddball.trialinfo == 3);

% Select Channel(s)
cfgs = [];
cfgs.channel = ROI;
roi = ft_selectdata(cfgs, clean_trials);

roi_erp = {};
n_trials = zeros([1 numel(cond_lab)]);
cfgavg = [];
cfgavg.keeptrials = 'yes';
cfgavg.trials = std_cdx;
roi_erp{1} = ft_timelockanalysis(cfgavg,roi);
ntrials(1) = size(roi_erp{1}.trial, 1);
cfgavg.trials = tar_cdx;
roi_erp{2} = ft_timelockanalysis(cfgavg,roi);
ntrials(2) = size(roi_erp{2}.trial, 1);
cfgavg.trials = odd_cdx;
roi_erp{3} = ft_timelockanalysis(cfgavg,roi);
ntrials(3) = size(roi_erp{3}.trial, 1);
