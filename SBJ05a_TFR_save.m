function SBJ05a_TFR_save(SBJ, proc_id, an_id)
% Filter data to create time-frequency representation (TFR):
%   Reconstruct and clean raw data, filter, cut trials to event, select channels, save
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
% OUTPUTS:
%   tfr [ft struct] - filtered data

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

% Load Behavioral Data
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

%% Load ICA, reconstruct and clean
if numel(SBJ_vars.block_name)>2; error('not ready for 3+ blocks!'); end
clean_data = fn_load_clean_experiment(SBJ, proc_id);

%% Determine Filter Padding
%   At a minimum, trial_lim_s must extend 1/2*max(filter_window) prior to 
%   the first data point to be estimated to avoid edges in the filtering window.
%   (ft_freqanalysis will return NaN for partially empty windows, e.g. an edge pre-trial,
%   but ft_preprocessing would return a filtered time series with an edge artifact.)
%   Also note this padding buffer should be at least 3x the slowest cycle
%   of interest.
if strcmp(cfg_tfr.method,'mtmconvol')
    pad_len = 0.5*max(cfg_tfr.t_ftimwin)*3;
elseif strcmp(cfg_tfr.method,'wavelet')
    % add 250 ms as a rule of thumb, or longer if necessary
    pad_len = 0.5*max([cfg_tfr.width/min(cfg_tfr.foi) 0.25]);
else
    error(['Unknown cfg_tfr.method: ' cfg_tfr.method]);
end
% Cut data to bsln_lim to be consistent across S and R locked (confirmed below)
%   Add extra 10 ms just because trimming back down to trial_lim_s exactly leave
%   one NaN on the end (smoothing that will NaN out everything)
trial_lim_s_pad = [min(an.bsln_lim)-pad_len an.trial_lim_s(2)+pad_len+0.01];

% Check that baseline will be included in data cut to trial_lim_s
if an.trial_lim_s(1) < an.bsln_lim(1)
    error(['ERROR: an.trial_lim_s does not include an.bsln_lim for an_id = ' an_id]);
end
% Check that trial_lim_s includes full baseline (e.g., zbtA)
if trial_lim_s_pad(2) < an.bsln_lim(2)+pad_len+0.01
    trial_lim_s_pad(2) = an.bsln_lim(2)+pad_len+0.01;
end

%% Obtain epoching events
% Load original trigger times for epoching
cfg_trl_unconcat = cell(size(SBJ_vars.block_name));
for b_ix = 1:numel(SBJ_vars.block_name)
    cfg = [];
    cfg.dataset             = SBJ_vars.dirs.raw_filename{b_ix};
    cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
    if strcmp(an.event_type,'S')
        cfg.trialdef.eventvalue = 1;%proc.event_code;
    elseif strcmp(an.event_type,'F')
        cfg.trialdef.eventvalue = 2;
    else
        error(['Unknown an.event_type: ' an.event_type]);
    end
    cfg.trialdef.prestim    = trial_lim_s_pad(1);
    cfg.trialdef.poststim   = trial_lim_s_pad(2);
    if startsWith(SBJ, 'EEG')
        cfg.tt_trigger_ix       = SBJ_vars.tt_trigger_ix;
        cfg.odd_trigger_ix      = SBJ_vars.odd_trigger_ix;
    end
    cfg.trialfun            = 'tt_trialfun';
    % Add downsample frequency since triggers are loaded from raw file
    cfg.resamp_freq         = proc.resample_freq;
    cfg_trl_unconcat{b_ix}  = ft_definetrial(cfg);
end

% Concatenate event times across blocks after adjusting times for gaps
cfg_trl = cfg_trl_unconcat{1};
if numel(SBJ_vars.block_name)>1
    % Get length of first block
    hdr = ft_read_header(SBJ_vars.dirs.raw_filename{1});
    endsample = hdr.nSamples;
    origFs = hdr.Fs;
    
    % Adjust event times for concatenated blocks
    cfg_trl_unconcat{2}.trl(:,1) = cfg_trl_unconcat{2}.trl(:,1)+endsample/origFs*proc.resample_freq;
    cfg_trl_unconcat{2}.trl(:,2) = cfg_trl_unconcat{2}.trl(:,2)+endsample/origFs*proc.resample_freq;
    cfg_trl.trl = vertcat(cfg_trl.trl, cfg_trl_unconcat{2}.trl);
end

% If the recording was started part way through, toss events not recorded
if any(cfg_trl.trl(:,1)<1)
    cfg_trl.trl(cfg_trl.trl(:,1)<1,:) = [];
end
event_onsets = cfg_trl.trl(:,1)-cfg_trl.trl(:,3);

% Cut trials
trials = ft_redefinetrial(cfg_trl, clean_data);

% Check trial_lim_s is within trial time (round to avoid annoying computer math)
if round(trial_lim_s_pad(1)+1/trials.fsample,3) < round(trials.time{1}(1),3) || ...
        round(trial_lim_s_pad(2)-1/trials.fsample,3) > round(trials.time{1}(end),3)
    warning('trial_lim_s_pad is outside data time bounds!');
end

%% Remove bad trials
% NOTE: exclude_trials from SBJ01 was not saved, so must be re-computed...
% Find trials that overlap with bad_epochs from raw visual inspection
load([SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat']);
if ~isempty(bad_epochs)
    bad_raw_trials = fn_find_trials_overlap_epochs(bad_epochs,1:size(clean_data.trial{1},2),...
        event_onsets,proc.trial_lim_s*clean_data.fsample);
else
    bad_raw_trials = [];
end

% Identify training and bad behavioral trials
%!!! error: the bhv fields don't match, so I somehow have old bhv files!
[bhv_orig] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);
training_ix = find(bhv_orig.blk==0);
rt_low_ix   = find(bhv_orig.rt <= proc.rt_bounds(1));
rt_high_ix  = find(bhv_orig.rt >= proc.rt_bounds(2));
exclude_trials = unique(vertcat(bad_raw_trials, training_ix, rt_low_ix, rt_high_ix));
fprintf(2,'\tWarning: Removing %i trials (%i bad_raw, %i training, %i rts)\n', numel(exclude_trials),...
    numel(bad_raw_trials), numel(training_ix), numel(rt_low_ix)+numel(rt_high_ix));

% Exclude original bad trials from SBJ02a (bad_epochs, training, RT outliers
%   Also select only channels in an.ROI!
cfgs = [];
cfgs.trials = setdiff([1:numel(trials.trial)], exclude_trials');
trials = ft_selectdata(cfgs, trials);

bhv_fields = fieldnames(bhv);
for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(exclude_trials) = [];
end

% Exclude variance-based outlier trials from SBJ02b/c
cfgs = [];
cfgs.trials  = setdiff([1:numel(clean_trials.trial)], SBJ_vars.trial_reject_ix);
cfgs.channel = an.ROI;
clean_trials = ft_selectdata(cfgs, clean_trials);

for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(SBJ_vars.trial_reject_ix) = [];
end

% Check that behavioral and EEG event triggers line up
if (numel(bhv.trl_n))~=numel(event_onsets)
    error(['Mismatch in behavioral and neural trial counts: ' num2str((numel(bhv.trl_n)))...
        ' behavioral; ' num2str(numel(event_onsets)) ' neural']);
end

%% Filter data
% all cfg_tfr options are specified in the an_vars
tfr_full = ft_freqanalysis(cfg_tfr, clean_data);

% Trim padding off
cfgs = [];
cfgs.latency = an.trial_lim_s;
tfr = ft_selectdata(cfgs, tfr_full);
clear tfr_full

%% Baseline Correction
switch an.bsln_type
    case {'zscore', 'demean', 'my_relchange'}
        tfr = fn_bsln_ft_tfr(tfr,an.bsln_lim,an.bsln_type,an.bsln_boots);
    case {'relchange','db'}
        cfgbsln = [];
        cfgbsln.baseline     = an.bsln_lim;
        cfgbsln.baselinetype = an.bsln_type;
        cfgbsln.parameter    = 'powspctrm';
        tfr = ft_freqbaseline(cfgbsln,tfr);
    otherwise
        error(['No baseline implemented for an.bsln_type: ' an.bsln_type]);
end

%% Save Results
data_out_fname = [SBJ_vars.dirs.proc SBJ '_' an_id '.mat'];
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'tfr');

end

