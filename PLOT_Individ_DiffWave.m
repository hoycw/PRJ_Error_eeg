function PLOT_Individ_DiffWave(SBJs, plt_id, an_id, fig_vis, save_fig, fig_ftype, conds)
%SBJS = cell of strings of subjects you want to plot {'EEG01', 'EEG02',
%etc}
%proc_id = 
%plt_id = 'ts_F15to28_evnts_sigPatch'
%an_id = 'ERP_Cz_F_trl15t28_flt05t20_stat06'
%fig_vis = whether or not the plot should pop up ('on' or 'off')
%save_fig = save  the figure or not (1 or 0)
%fig_ftype = file type the figure should be saved as (i.e. 'png')
%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);
addpath(ft_dir);
ft_defaults
%% Load preprocessed data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat'];
load(clean_bhv_fname);
data_cleanname = [SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat'];
load(data_cleanname);
%% Compute ERPs
% Select Channel(s)
cfgs = [];
cfgs.channel = an.ROI;
roi = ft_selectdata(cfgs, clean_trials);

cond_lab = {'easy', 'hard'};
cond_lab_hit = [1 0];
cond_colors = {[0.6350, 0.0780, 0.1840], [0.3010, 0.7450, 0.9330], 	[0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560]}
num_conds = (size(cond_lab, 2)+size(cond_lab_hit, 2));
roi_erp  = cell(num_conds);
n_trials = zeros(1, num_conds);
cfgavg = [];
cfgavg.keeptrials = 'yes';

