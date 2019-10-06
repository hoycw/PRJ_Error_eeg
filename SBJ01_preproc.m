function SBJ01_preproc(SBJ, proc_id)
%% Import and preprocess EEG data. Null the bad channels and null the marked bad epochs from SBJ00.  Cut the data into trials and then run ICA.
% INPUTS:
%   SBJ [str] - name of the subject to load
%   proc_id [str] - name of processing pipeline

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Set up processing and SBJ variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' odd_proc_id '_vars.m'];
eval(proc_vars_cmd);

%% Load and preprocess the data
% Update ear ref channels with prefix/suffix
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

% Load and preprocess
cfg = [];
cfg.continuous = 'yes'; 
% cfg.lpfilter   = proc.lp_yn;
% cfg.hpfilter   = proc.hp_yn;
cfg.bpfilter   = proc.bp_yn;
cfg.bpfreq     = proc.bp_freq;
cfg.bpfiltord  = proc.bp_order;
cfg.demean     = proc.demean_yn;
cfg.reref      = proc.reref_yn;
cfg.refmethod  = proc.ref_method;
cfg.refchannel = {ear_lab1, ear_lab2};

data = cell(size(SBJ_vars.block_name));
for b_ix = 1:numel(SBJ_vars.block_name)
    cfg.dataset    = SBJ_vars.dirs.raw_filename{b_ix};
    data{b_ix} = ft_preprocessing(cfg);
end
data = fn_concat_blocks(data);

%% Downsample
if strcmp(proc.resample_yn,'yes')
    cfg = [];
    cfg.resamplefs = proc.resample_freq;
    data = ft_resampledata(cfg, data);
end

%% Fix channel labels
% Remove bad channels
null_neg = cell(size(SBJ_vars.ch_lab.null));
for null_ix = 1:numel(SBJ_vars.ch_lab.null)
    null_neg{null_ix} = ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.null{null_ix} SBJ_vars.ch_lab.suffix];
end
bad_neg = cell(size(SBJ_vars.ch_lab.bad));
for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
    bad_neg{bad_ix} = ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix];
end
rep_neg = cell(size(SBJ_vars.ch_lab.replace));
for rep_ix = 1:numel(SBJ_vars.ch_lab.replace)
    rep_neg{rep_ix} = ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.replace{rep_ix}{1} SBJ_vars.ch_lab.suffix];
end
ears_neg = {['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{1} SBJ_vars.ch_lab.suffix],...
            ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{2} SBJ_vars.ch_lab.suffix]};

cfg = [];
cfg.channel = [{'all'}, bad_neg, rep_neg, ears_neg, null_neg];
data = ft_selectdata(cfg,data);

% Strip Pre/Suffix if Necessary
for ch_ix = 1:numel(data.label)
    data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
    data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
    % Replace label of the externals with actual electrode name
    for x = 1:numel(SBJ_vars.ch_lab.replace)
        if strcmp(SBJ_vars.ch_lab.replace{x}{2},data.label{ch_ix})
            data.label{ch_ix} = SBJ_vars.ch_lab.replace{x}{1};
        end
    end
end

%% Extract and process EOG
% Grab first horizontal channel
cfg = [];
cfg.channel = SBJ_vars.ch_lab.eog_h{1};
eog_h = ft_selectdata(cfg,data);
eog_h.label{1} = 'eog_h';
% Subtract second horizontal channel
eog_h2_ix = strcmp(data.label,SBJ_vars.ch_lab.eog_h{2});
eog_h.trial{1}(1,:) = eog_h.trial{1}(1,:)-data.trial{1}(eog_h2_ix,:);

% Grab first vertical channel
cfg.channel = SBJ_vars.ch_lab.eog_v{1};
eog_v = ft_selectdata(cfg,data);
eog_v.label{1} = 'eog_v';
% Subtract second vertical channel
eog_v2_ix = strcmp(data.label,SBJ_vars.ch_lab.eog_v{2});
eog_v.trial{1}(1,:) = eog_v.trial{1}(1,:)-data.trial{1}(eog_v2_ix,:);

% Combine bipolar EOG
cfg = [];
eog = ft_appenddata(cfg,eog_h,eog_v);

% Remove unipolar EOG
warning('WARNING!!! Assuming Fp2 is the second vertical EOG!');
eog_v_low_ix = ~strcmp(SBJ_vars.ch_lab.eog_v,'Fp2');    % only toss the lower one
eog_neg = [fn_ch_lab_negate(SBJ_vars.ch_lab.eog_h),fn_ch_lab_negate(SBJ_vars.ch_lab.eog_v(eog_v_low_ix))];
trig_neg = fn_ch_lab_negate({SBJ_vars.ch_lab.trigger});
cfg = [];
cfg.channel = [{'all'}, eog_neg, trig_neg];
data = ft_selectdata(cfg,data);

%% Cut into trials
% % Must segment before downsampling because trigger channel read from
% % original file
% % Need to cut into trials here so that can NaN out bad trials later in
% % script
% cfg = [];
% cfg.dataset             = SBJ_vars.dirs.raw_filename;
% cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
% cfg.trialdef.eventvalue = proc.event_code;        % feedback cocde
% cfg.trialdef.prestim    = proc.trial_lim_s(1);
% cfg.trialdef.poststim   = proc.trial_lim_s(2);
% cfg.trialfun            = 'tt_trialfun';%'ft_trialfun_general';%
% cfg_trl = ft_definetrial(cfg);
% 
% % If the recording was started part way through, toss events not recorded
% if any(cfg_trl.trl(:,1)<1)  
%     cfg_trl.trl(cfg_trl.trl(:,1)<1,:) = [];
% end
% event_onsets = cfg_trl.trl(:,1)-cfg_trl.trl(:,3);
% 
% trials = ft_redefinetrial(cfg_trl,data);
% eog_trials = ft_redefinetrial(cfg_trl,eog);


%% ICA
% Load raw bad epochs to NaN out
bad_epochs = fn_combine_raw_bad_epochs(SBJ);

% Preprocess data for ICA (NaN bad epochs)
if ~isempty(bad_epochs)
    cfg = [];
    cfg.artfctdef.visual.artifact = bad_epochs;
    cfg.artfctdef.reject          = 'nan';
    data = ft_rejectartifact(cfg, data);
end
if strcmp(proc.ICA_hp_yn,'yes')
    cfg = [];
    cfg.hpfilter = proc.ICA_hp_yn;
    cfg.hpfreq   = proc.ICA_hp_freq;
    data_ICA = ft_preprocessing(cfg, data);
else
    data_ICA = data;
end

% Run ICA
cfg = [];
cfg.method = 'runica'; % default, uses the implementation from EEGLAB
icomp = ft_componentanalysis(cfg, data_ICA);

%% Save ICA unmixing and topolabel, data, eog, and bad_epochs
icaunmixing  = icomp.unmixing;
icatopolabel = icomp.topolabel;

data_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
save(data_fname, 'icaunmixing', 'icatopolabel', 'data', 'eog','bad_epochs');


