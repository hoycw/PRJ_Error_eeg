function TFR_preproc(SBJ, proc_id)
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/', ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);
addpath(ft_dir);
ft_defaults

%% Processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

%% Load data
% Load EEG
load([SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat']);
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat']);
% Load Behavior
[bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);

%% Remove time range of bad trials
%% Get times of trials
% Need to recut trials on updated data with the nans
for b_ix = 1:numel(SBJ_vars.block_name)
cfg = [];
cfg.dataset             = SBJ_vars.dirs.raw_filename{b_ix};
cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
cfg.trialdef.eventvalue = proc.event_code;        % feedback cocde
cfg.trialdef.prestim    = proc.trial_lim_s(1);
cfg.trialdef.poststim   = proc.trial_lim_s(2);
if startsWith(SBJ, 'EEG')
cfg.tt_trigger_ix       = SBJ_vars.tt_trigger_ix;
cfg.odd_trigger_ix      = SBJ_vars.odd_trigger_ix;
end
cfg.trialfun            = 'tt_trialfun';
% Add downsample frequency since triggers are loaded from raw file
cfg.resamp_freq         = proc.resample_freq;
cfg_trl_unconcat{b_ix}  = ft_definetrial(cfg);
end
hdr = ft_read_header(SBJ_vars.dirs.raw_filename{1});
endsample = hdr.nSamples;
origFs = hdr.Fs;
if numel(SBJ_vars.block_name)>1
for b_ix = 2:numel(SBJ_vars.block_name)
cfg_trl_unconcat{b_ix}.trl(:,1) = cfg_trl_unconcat{b_ix}.trl(:,1)+endsample/origFs*proc.resample_freq;
cfg_trl_unconcat{b_ix}.trl(:,2) = cfg_trl_unconcat{b_ix}.trl(:,2)+endsample/origFs*proc.resample_freq;
cfg_trl_unconcat{1}.trl = vertcat(cfg_trl_unconcat{b_ix-1}.trl, cfg_trl_unconcat{b_ix}.trl);
end
end
cfg_trl = cfg_trl_unconcat{1};
% If the recording was started part way through, toss events not recorded
if any(cfg_trl.trl(:,1)<1)
cfg_trl.trl(cfg_trl.trl(:,1)<1,:) = [];
end
event_onsets = cfg_trl.trl(:,1)-cfg_trl.trl(:,3);

if ~isempty(bad_epochs)
bad_raw_trials = fn_find_trials_overlap_epochs(bad_epochs,1:size(data.trial{1},2),...
                      event_onsets,proc.trial_lim_s*data.fsample);
else
bad_raw_trials = [];
end
bad_times = cell(numel(bad_raw_trials));
bad_trial_times = zeros(numel(bad_raw_trials), 2);
good_trial_times = [1: data.sampleinfo(2)];
for trial_ix = 1: numel(bad_raw_trials)
    bad_trial_times(trial_ix,1) = cfg_trl.trl(bad_raw_trials(trial_ix), 1);
    bad_trial_times(trial_ix,2) = cfg_trl.trl(bad_raw_trials(trial_ix), 2);
    bad_times{trial_ix} = [bad_trial_times(trial_ix,1):bad_trial_times(trial_ix,2)];
    good_trial_times = setdiff(good_trial_times, bad_times{trial_ix});
end
cfg = [];
cfg.trials = good_trial_times; %How to do this?
%data = ft_selectdata(cfg, data);
% Rebuild the components
cfg           = [];
cfg.unmixing  = icaunmixing;
cfg.topolabel = icatopolabel;
ica           = ft_componentanalysis(cfg, data);
%% IC rejection
cfg = [];
cfg.component = unique([SBJ_vars.ica_reject, heog_ics, veog_ics]);
cfg.demean = 'no';
clean_trials = ft_rejectcomponent(cfg, ica);

%% Repair Bad Channels
%Adding them back in enables ft_databrowser to plot full cap correctly
cfg = [];
cfg.method         = 'average';
cfg.missingchannel = SBJ_vars.ch_lab.bad(:); % not in data (excluded from ica)
cfg.layout         = 'biosemi64.lay';

cfgn = [];
cfgn.channel = 'all';
cfgn.layout  = 'biosemi64.lay';
cfgn.method  = 'template';
cfg.neighbours = ft_prepare_neighbours(cfgn);

clean_trials = ft_channelrepair(cfg, clean_trials);
%% Filter by TFR

%% Cut into trials -- FIX THIS: I don't think will work if you already cut out the trials!
% Need to recut trials on updated data with the nans
for b_ix = 1:numel(SBJ_vars.block_name)
cfg = [];
cfg.dataset             = SBJ_vars.dirs.raw_filename{b_ix};
cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
cfg.trialdef.eventvalue = proc.event_code;        % feedback cocde
cfg.trialdef.prestim    = proc.trial_lim_s(1);
cfg.trialdef.poststim   = proc.trial_lim_s(2);
if startsWith(SBJ, 'EEG')
cfg.tt_trigger_ix       = SBJ_vars.tt_trigger_ix;
cfg.odd_trigger_ix      = SBJ_vars.odd_trigger_ix;
end
cfg.trialfun            = 'tt_trialfun';
% Add downsample frequency since triggers are loaded from raw file
cfg.resamp_freq         = proc.resample_freq;
cfg_trl_unconcat{b_ix}  = ft_definetrial(cfg);
end
hdr = ft_read_header(SBJ_vars.dirs.raw_filename{1});
endsample = hdr.nSamples;
origFs = hdr.Fs;
if numel(SBJ_vars.block_name)>1
for b_ix = 2:numel(SBJ_vars.block_name)
cfg_trl_unconcat{b_ix}.trl(:,1) = cfg_trl_unconcat{b_ix}.trl(:,1)+endsample/origFs*proc.resample_freq;
cfg_trl_unconcat{b_ix}.trl(:,2) = cfg_trl_unconcat{b_ix}.trl(:,2)+endsample/origFs*proc.resample_freq;
cfg_trl_unconcat{1}.trl = vertcat(cfg_trl_unconcat{b_ix-1}.trl, cfg_trl_unconcat{b_ix}.trl);
end
end
cfg_trl = cfg_trl_unconcat{1};
% If the recording was started part way through, toss events not recorded
if any(cfg_trl.trl(:,1)<1)
cfg_trl.trl(cfg_trl.trl(:,1)<1,:) = [];
end
event_onsets = cfg_trl.trl(:,1)-cfg_trl.trl(:,3);

% Cut the data into trials
trials = ft_redefinetrial(cfg_trl,data);
eog_trials = ft_redefinetrial(cfg_trl,eog);

% Check that behavioral and EEG event triggers line up
if (numel(bhv.trl_n))~=numel(event_onsets)
error(['Mismatch in behavioral and neural trial counts: ' num2str((numel(bhv.trl_n)))...
' behavioral; ' num2str(numel(event_onsets)) ' neural']);
end
end
