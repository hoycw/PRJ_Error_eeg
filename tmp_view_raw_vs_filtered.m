if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Load raw (from EEG00)
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

cfg=[];
cfg.dataset  = SBJ_vars.dirs.raw_filename;
% cfg.demean   = 'yes';
% cfg.lpfilter = 'no';
cfg.hpfilter = 'yes';
cfg.hpfreq   = 1; %changed to 1 because EP07 was erroring line 54 of filter_with_correlation 
% cfg.bpfilter = 'yes';
% cfg.bpfreq   = [0.5 40];%0.1 is too low for filter settings, 20 is too low to see muscle artifact, consider ditching filtering?
raw = ft_preprocessing(cfg);

% Load cfg with plotting parameters
load([root_dir 'PRJ_Error_eeg/scripts/utils/cfg_plot_eeg.mat']);

% Load previous if available
out_fname = [SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat'];
load(out_fname);
cfg_plot.artfctdef.visual.artifact = bad_epochs;

if isfield(SBJ_vars.ch_lab,'prefix')
    ear_lab1 = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{1}];
    ear_lab2 = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{2}];
else
    ear_lab1 = [SBJ_vars.ch_lab.ears{1}];
    ear_lab2 = [SBJ_vars.ch_lab.ears{2}];
end
if isfield(SBJ_vars.ch_lab,'suffix')
    ear_lab1 = [ear_lab1 SBJ_vars.ch_lab.suffix];
    ear_lab2 = [ear_lab2 SBJ_vars.ch_lab.suffix];
end
bad_neg = {};
for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
    bad_neg = {bad_neg{:},['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix]};
end
rep_neg = {};
for rep_ix = 1:numel(SBJ_vars.ch_lab.replace)
    rep_neg = {rep_neg{:},['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.replace{rep_ix}{1} SBJ_vars.ch_lab.suffix]};
end
ears_neg = {['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{1} SBJ_vars.ch_lab.suffix],...
            ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{2} SBJ_vars.ch_lab.suffix]};
cfg = [];
cfg.channel = {'all',bad_neg{:},rep_neg{:},ears_neg{:},'-Status'};

raw = ft_selectdata(cfg,raw);

%% Load preproc (from EEG02)
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);

data_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
load(data_fname);

%% plot
ft_databrowser(cfg_plot,raw);
ft_databrowser(cfg_plot,data);

%% stupdi bad trials line up
% times = [1971 1981];
% trial_ix = zeros(size(times));
% for t = 1:numel(times)
%     tmp = find(raw.time{1}>times(t));
%     tmp2 = find(clean_trials.sampleinfo(:,1)>tmp(1));
%     trial_ix(t) = tmp2(1);
%     fprintf('time = %.02f; trial_ix = %d = %.02f - %.02f\n',times(t),trial_ix(t),clean_trials.sampleinfo(trial_ix(t),:)/clean_trials.fsample);
% end
trial_ix = zeros([size(bad_epochs,1) 1]);
for t = 1:size(bad_epochs,1)
%     tmp = find(raw.time{1}>times(t));
    tmp2 = find(clean_trials.sampleinfo(:,1)>bad_epochs(t,1));
    trial_ix(t) = tmp2(1);
    fprintf('time = %.02f - %.02f; trial_ix = %d = %.02f - %.02f\n',bad_epochs(t,:)/raw.fsample,trial_ix(t),clean_trials.sampleinfo(trial_ix(t),:)/clean_trials.fsample);
end

cfgp = [];
cfgp.trial = trial_ix;
cfgp.viewmode = 'vertical';
% bad_trials = ft_