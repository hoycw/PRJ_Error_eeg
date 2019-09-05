function EEG02c_trial_rejection(SBJ,proc_id, visual)
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
data_cleanname = [SBJ_vars.dirs.preproc SBJ '_clean02b_' proc_id '.mat'];
load(data_cleanname)
clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav02b_' proc_id '_clean.mat'];
load(clean_bhv_fname);

%% Eliminate trials
    cfgs = [];
    trialsnum_oddball = numel(bhv.trl_n) - bhv.numtrials;
    bhv.numtrials = bhv.numtrials - numel(find(SBJ_vars.trial_reject_ix>bhv.numtrials_odd));
    bhv.numtrials_odd = bhv.numtrials_odd - numel(find(SBJ_vars.trial_reject_ix<=bhv.numtrials_odd));
    cfgs.trials = setdiff([1:numel(clean_trials.trial)], SBJ_vars.trial_reject_ix); % note: trial_reject_ix is the index of the values in clean trials that show up in data_browser.
    %WARNING!! If you do databrowser after doing the steps above and
    %deciding to reject trials in the gooey of lines 44 - 52, your indices
    %will be incorrect.  they need to be the indices of clean_trials in line
    %40. 
    clean_trials = ft_selectdata(cfgs, clean_trials);
%% FINAL CHECK
% Load cfg with plotting parameters
if visual
    cfg_plot.ylim = [-15 15];  %this is different from initial but its in order to see anything unusual
    cfg_plot.viewmode = 'vertical';
    ft_databrowser(cfg_plot, clean_trials);
end

for f_ix = 1:numel(bhv_fields)
    if ~(strcmp(bhv_fields{f_ix}, 'numtrials')) && ~(strcmp(bhv_fields{f_ix}, 'numtrials_odd'))
        bhv.(bhv_fields{f_ix})(SBJ_vars.trial_reject_ix) = [];
    end
end
%Break the trials up

cfg_oddball = [];
totalnum = numel(bhv.trl_n);
cfg_oddball.trials = [1:bhv.numtrials_odd];
clean_oddball = ft_selectdata(cfg_oddball, clean_trials);
cfg_trials.trials = [bhv.numtrials_odd + 1:totalnum];
clean_trials = ft_selectdata(cfg_trials, clean_trials);

a = find(bhv.trl_n == 1);
real_trial_start = max(a)-1;
for f_ix = 1:numel(bhv_fields)
    if ~(strcmp(bhv_fields{f_ix}, 'numtrials')) && ~(strcmp(bhv_fields{f_ix}, 'numtrials_odd'))
        bhv_oddball.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(1:real_trial_start);
    end
end



%% Save outputs
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat'];
save(clean_data_fname, '-v7.3', 'clean_trials');

clean_data_oddball_fname = [SBJ_vars.dirs.preproc SBJ '_clean_oddball_' proc_id '.mat'];
save(clean_data_oddball_fname, '-v7.3', 'clean_oddball');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv');

clean_bhv_oddball_fname = [SBJ_vars.dirs.events SBJ '_behav_oddball_' proc_id '_clean.mat'];
save(clean_bhv_oddball_fname, '-v7.3', 'bhv_oddball');

