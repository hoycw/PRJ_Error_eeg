function SBJ09a_TFR_stats(SBJ,conditions,proc_id,an_id)
error('switch to save styel then run separate stats scripts like SBJ08a!');
% Calculates time frequency representation, computes cluster-based statistics, and plots the results
% INPUTS:
%   SBJ [str] - dataset to be processed
%   conditions [str] - experimental conditions to be compared
%   proc_id [str] - which processed data pipeline to get the data
%   an_id [str] - which analysis variables to use

% clear all; %close all;
%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
% Load Data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

% Select Conditions of Interest
[cond_lab, ~, ~] = fn_condition_label_styles(conditions);
cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n));
end

%% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Cut into Trials
% Pad trial_lim_s by 1/2 lowest frequency window length to avoid NaNs in epoch of interest
% Add 10 ms just because trimming back down to trial_lim_s exactly leaves
% one NaN on the end, so smoothing will NaN out everything

%!!! BEWARE: This is suboptimal for R-locked because trial_lim_s(1)=-0.5,
%which adds 250ms to the time series that isn't necessary; the realign_tfr
%function should still cut to the desired data, but it'll take longer.
if strcmp(cfg_tfr.method,'mtmconvol')
    trial_lim_s_pad = [an.trial_lim_s(1)-max(cfg_tfr.t_ftimwin)/2 ...
                       an.trial_lim_s(2)+max(cfg_tfr.t_ftimwin)/2+0.01];
elseif strcmp(cfg_tfr.method,'wavelet')
    trial_lim_s_pad = [an.trial_lim_s(1)-cfg_tfr.width/min(an.foi_center)*2 ...
                       an.trial_lim_s(2)+cfg_tfr.width/min(an.foi_center)*2+0.01];
else
    error(['Cannot adjust baseline time window, unknown TFR filtering method: ' cfg_tfr.method]);
end

% Always normalize to pre-stimulus baseline for HFA
bsln_events = trial_info.word_onset;
if strcmp(an.evnt_lab,'S')
    % Cut to desired trial_lim_s
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,trial_info.condition_n',...
        round(trial_lim_s_pad*roi.fsample));
elseif strcmp(an.evnt_lab,'R')
    % Check that baseline will be included in trial_lim_s
    if an.trial_lim_s(1)>an.bsln_lim(1)
        error(['ERROR: an.trial_lim_s does not include an.bsln_lim for an_id = ' an_id]);
    end
    % Cut out to max_RT+trial_lim_s(2)+max(cfg_hfa.t_ftimwin)
    max_RT  = max(trial_info.response_time);
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) max_RT+trial_lim_s_pad(2)]*roi.fsample));
else
    error(['Unknown evnt_lab: ' an.evnt_lab]);
end

%% Compute TFRs
% all cfg_tfr options are specified in the an_vars
tfr      = {};
n_trials = zeros([1 numel(cond_lab)]);
for cond_ix = 1:numel(cond_lab)
    fprintf('===================================================\n');
    fprintf('------------- TFR Calculations for %s ----------\n',cond_lab{cond_ix});
    fprintf('===================================================\n');
    cfg_tfr.trials = find(cond_idx(cond_ix,:)==1);
    tfr{cond_ix}   = ft_freqanalysis(cfg_tfr, roi_trl);
    
    % Trim back down to original trial_lim_s to exclude NaNs
    if strcmp(an.evnt_lab,'S')
        cfg_trim = [];
        cfg_trim.latency = an.trial_lim_s;
        tfr{cond_ix} = ft_selectdata(cfg_trim,tfr{cond_ix});
    elseif strcmp(an.evnt_lab,'R')
        cfg_trim = [];
        cfg_trim.latency = [an.trial_lim_s(1) max_RT+an.trial_lim_s(2)];
        tfr{cond_ix} = ft_selectdata(cfg_trim,tfr{cond_ix});
    else
        error(['Unknown evnt_lab: ' an.evnt_lab]);
    end
    
    % Grab n_trials for design matrix
    n_trials(cond_ix) = size(tfr{cond_ix}.trialinfo,1);
end

%% Baseline Correction
fprintf('===================================================\n');
fprintf('------------- Baseline Correction for %s ----------\n',cond_lab{cond_ix});
fprintf('===================================================\n');
for cond_ix = 1:numel(cond_lab)
    switch an.bsln_type
        case {'zscore', 'demean', 'my_relchange'}
            tfr{cond_ix} = fn_bsln_ft_tfr(tfr{cond_ix},an.bsln_lim,an.bsln_type,an.bsln_boots);
        case 'relchange'
            cfgbsln = [];
            cfgbsln.baseline     = an.bsln_lim;
            cfgbsln.baselinetype = an.bsln_type;
            cfgbsln.parameter    = 'powspctrm';
            tfr{cond_ix} = ft_freqbaseline(cfgbsln,tfr{cond_ix});
        otherwise
            error(['No baseline implemented for an.bsln_type: ' an.bsln_type]);
    end
end

%% Re-align to event of interest if necessary (e.g., response)
if strcmp(an.evnt_lab,'R')
    for cond_ix = 1:numel(cond_lab)
        tfr{cond_ix} = fn_realign_tfr_s2r(tfr{cond_ix},trial_info.response_time,an.trial_lim_s);
    end
elseif ~strcmp(an.evnt_lab,'S')
    error(['ERROR: unknown evnt_lab ' an.evnt_lab]);
end

%% Run Statistics
fprintf('===================================================\n');
fprintf('--------------------- Statistics ------------------\n');
fprintf('===================================================\n');
% Create design matrix
design = zeros(2,sum(n_trials));
for an_ix = 1:numel(cond_lab)
    if an_ix==1
        design(1,1:n_trials(an_ix)) = an_ix;                                % Conditions (Independent Variable)
        design(2,1:n_trials(an_ix)) = 1:n_trials(an_ix);                    % Trial Numbers
    else
        design(1,sum(n_trials(1:an_ix-1))+1:sum(n_trials(1:an_ix)))= an_ix; % Conditions (Independent Variable)
        design(2,sum(n_trials(1:an_ix-1))+1:sum(n_trials(1:an_ix)))= 1:n_trials(an_ix);
    end
end

% Calculate statistics
cfg_stat.design           = design;
[stat] = ft_freqstatistics(cfg_stat, tfr{:});

%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',data_out_filename);
fprintf('===================================================\n');
save(data_out_filename, '-v7.3','stat','tfr','an');

end
