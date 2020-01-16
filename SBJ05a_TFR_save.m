function SBJ05a_TFR_save(SBJ, proc_id, an_id)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Load Data
%!!! load SBJ01 uncut output
%load([SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

%% Re-clean uncut data
%!!! Sheila function to project to ICA space, bring it back without bad
%components
clean_data = fn_TFR_clean(SBJ, proc_id);
%% Re-cut into trials
%!!!Colin will decide whether to write this new or modify realignment below
% Sheila: I left it so that the bad trials get rejected here
%% Realign data to desired event
if ~strcmp(proc.event_type,an.event_type)
    cfg = [];
    % Match desired time to closest sample index
    if strcmp(proc.event_type,'S') && strcmp(an.event_type,'F')
        prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);
        cfg.offset = -(prdm_vars.target + prdm_vars.fb_delay)*clean_trials.fsample;
    elseif strcmp(proc.event_type,'S') && strcmp(an.event_type,'R')
        cfg.offset = -bhv.rt*clean_trials.fsample;
    elseif strcmp(proc.event_type,'F')
        error('F-locked preprocessing can only be used for F-locked analysis!');
    elseif strcmp(proc.event_type,'R')% && strcmp(an.event_type,'S')
        error('Why were you doing R-locked preprocessing?');
        %error('cannot do S-locked analysis with R-locked data!');
    else
        error('unknown combination of proc and an event_types');
    end
    % Convert time axis to new event:
    %   basically: data.time{i} = data.time{i} + offset(i)/data.fsample;
    %   therefore, negative offset will shift time axis "back"
    roi = ft_redefinetrial(cfg, clean_trials);
else
    roi = clean_trials;
end

%% Select Data for TFR
% Pad trial_lim_s by 1/2 lowest frequency window length to avoid NaNs in epoch of interest
% Add 10 ms just because trimming back down to trial_lim_s exactly leaves
% one NaN on the end, so smoothing will NaN out everything
if strcmp(cfg_tfr.method,'mtmconvol')
    trial_lim_s_pad = [an.trial_lim_s(1)-max(cfg_tfr.t_ftimwin)/2 ...
                       an.trial_lim_s(2)+max(cfg_tfr.t_ftimwin)/2+0.01];
elseif strcmp(cfg_tfr.method,'wavelet')
    trial_lim_s_pad = [an.trial_lim_s(1)-cfg_tfr.width/min(cfg_tfr.foi)*2 ...
                       an.trial_lim_s(2)+cfg_tfr.width/min(cfg_tfr.foi)*2+0.01];
else
    error(['Cannot adjust baseline time window, unknown TFR filtering method: ' cfg_tfr.method]);
end

% Check window consistency
% Check trial_lim_s is within trial time (round to avoid annoying computer math)
if round(trial_lim_s_pad(1)+1/roi.fsample,3) < round(roi.time{1}(1),3) || ...
        round(trial_lim_s_pad(2)-1/roi.fsample,3) > round(roi.time{1}(end),3)
    warning('trial_lim_s_pad is outside data time bounds!');
end

%% Select Channel(s)
cfgs = [];
cfgs.channel = an.ROI;
cfgs.latency = trial_lim_s_pad;
roi = ft_selectdata(cfgs, roi);

%% Compute TFRs
% all cfg_tfr options are specified in the an_vars
tfr  = ft_freqanalysis(cfg_tfr, roi);

% Trim padding off
cfgs = [];
cfgs.latency = an.trial_lim_s;
tfr = ft_selectdata(cfgs, tfr);

%% Baseline Correction
switch an.bsln_type
    case {'zscore', 'demean', 'my_relchange'}
        tfr = fn_bsln_ft_tfr(tfr,an.bsln_lim,an.bsln_type,an.bsln_boots);
    case 'relchange'
        cfgbsln = [];
        cfgbsln.baseline     = an.bsln_lim;
        cfgbsln.baselinetype = an.bsln_type;
        cfgbsln.parameter    = 'powspctrm';
        tfr = ft_freqbaseline(cfgbsln,tfr);
    otherwise
        error(['No baseline implemented for an.bsln_type: ' an.bsln_type]);
end

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'tfr');
end

