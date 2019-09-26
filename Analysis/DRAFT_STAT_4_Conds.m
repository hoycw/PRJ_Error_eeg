function DRAFT_STAT_4_Conds(SBJs, plt_id, an_id, fig_vis, save_fig, fig_ftype)
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
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);
addpath(ft_dir);
ft_defaults
%% Load preprocessed data
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cond_lab_hit = [1 0];
cond_lab = {'easy', 'hard'};
cond_colors = {[0.6350, 0.0780, 0.1840], [0.3010, 0.7450, 0.9330], 	[0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560]};
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/Grps/' an_id '_4_Conds/'];
clean_stats_fname = [fig_dir an_id 'GRP_4_Conds.mat'];
load(clean_stats_fname);
%%Find peak:
[~, gstat_bounds(1)] = min(abs(gstat{2}.time-1.8));
[~, gstat_bounds(2)] = min(abs(gstat{2}.time-2.07));
[peak, index] = findpeaks(gstat{2}.avg(gstat_bounds(1):gstat_bounds(2)))
index_pos = min(index);



