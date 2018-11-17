function EEG01_preproc(SBJ, proc_id)
%% Import and preprocess EEG data, then run ICA

%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath(genpath([root_dir 'PRJ_Error_eeg/scripts/']));
addpath(ft_dir);
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);

%% Load, preprocess, and downsample the data
% Update ears ref channels with prefix/suffix
if isfield(SBJ_vars.ch_lab,'prefix')
    ear_lab1 = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{1}];
    ear_lab2 = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{2}];
else
    ear_lab1 = [SBJ_vars.ch_lab.ears{1}];
    ear_lab1 = [SBJ_vars.ch_lab.ears{2}];
end
if isfield(SBJ_vars.ch_lab,'suffix')
    ear_lab1 = [ear_lab1 SBJ_vars.ch_lab.suffix];
    ear_lab2 = [ear_lab2 SBJ_vars.ch_lab.suffix];
end

% Load and preprocess
cfg=[];
cfg.dataset    = SBJ_vars.dirs.raw_filename;
cfg.continuous = 'yes'; %!!! try segmenting trial on import
cfg.lpfilter   = proc_vars.lp_yn;
cfg.hpfilter   = proc_vars.hp_yn;
cfg.bpfilter   = proc_vars.bp_yn;
cfg.bpfreq     = proc_vars.bp_freq;
cfg.demean     = proc_vars.demean_yn;
cfg.reref      = proc_vars.reref_yn;
cfg.refmethod  = proc_vars.ref_method;
cfg.refchannel = {ear_lab1, ear_lab2};
data = ft_preprocessing(cfg);

%% Fix channel labels
% Remove bad channels
bad_neg = {};
for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
    bad_neg = {bad_neg{:},['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix]};
end
rep_neg = {};
for rep_ix = 1:numel(SBJ_vars.ch_lab.replacements)
    rep_neg = {rep_neg{:},['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.replacements{rep_ix}{1} SBJ_vars.ch_lab.suffix]};
end
ears_neg = {['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{1} SBJ_vars.ch_lab.suffix],...
            ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{2} SBJ_vars.ch_lab.suffix]};
cfg = [];
cfg.channel = {'all',bad_neg{:},rep_neg{:},ears_neg{:}};
data = ft_selectdata(cfg,data);

% Strip Pre/Suffix if Necessary
for ch_ix = 1:numel(data.label)
    data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
    data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
    for x = 1:numel(SBJ_vars.ch_lab.replacements)
        if strcmp(SBJ_vars.ch_lab.replacements{x}{2},data.label{ch_ix})
            data.label{ch_ix} = SBJ_vars.ch_lab.replacements{x}{1}; % replaces the label of the externals with the channel they represent
        end
    end
end

%% Downsample
if strcmp(proc_vars.resample_yn,'yes')
    cfg=[];
    cfg.resamplefs = proc_vars.resample_freq;
    data = ft_resampledata(cfg, data);
end

%% Extract and process EOG
% Grab first channel
cfg = [];
cfg.channel = SBJ_vars.ch_lab.eog_h{1};
eog_h = ft_selectdata(cfg,data);
eog_h.label{1} = 'eog_h';
% Subtract second channel
eog_h2_ix = strcmp(data.label,SBJ_vars.ch_lab.eog_h{2});
eog_h.trial{1}(1,:) = eog_h.trial{1}(1,:)-data.trial{1}(eog_h2_ix,:);

% Grab first channel
cfg.channel = SBJ_vars.ch_lab.eog_v{1};
eog_v = ft_selectdata(cfg,data);
eog_v.label{1} = 'eog_v';
% Subtract second channel
eog_v2_ix = strcmp(data.label,SBJ_vars.ch_lab.eog_v{2});
eog_v.trial{1}(1,:) = eog_v.trial{1}(1,:)-data.trial{1}(eog_v2_ix,:);

% Combine bipolar EOG
cfg = [];
eog = ft_appenddata(cfg,eog_h,eog_v);

% Remove unipolar EOG
eog_v_low_ix = ~strcmp(SBJ_vars.ch_lab.eog_v,'Fp2');    % only toss the lower one
eog_neg = [fn_ch_lab_negate(SBJ_vars.ch_lab.eog_h),fn_ch_lab_negate(SBJ_vars.ch_lab.eog_v(eog_v_low_ix))];
cfg = [];
cfg.channel = {'all',eog_neg{:}};
data = ft_selectdata(cfg,data);

% Baseline-correction options- Not doing right now
%cfg.demean          = 'yes'; %baseline correct data before you resample?
%cfg.baselinewindow  = [-0.25 -0.05]; %set baseline to -.2 to 0 seconds, then correct it before you resample

%% Cut into trials
cfg = [];
cfg.dataset             = SBJ_vars.dirs.raw_filename;
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = 2;        % feedback cocde
cfg.trialdef.prestim    = -0.7;
cfg.trialdef.poststim   = 1.5;
cfg.trialfun            = 'ft_trialfun_general';
cfg = ft_definetrial(cfg);
trials = ft_preprocessing(cfg,data);
eog_trials = ft_preprocessing(cfg,eog);

%!!! WTF was this about?
% cfg = [];
% cfg.latency = [-0.2 1];
% trials = ft_selectdata(cfg,trials);

% cfg.artfctdef.visual.artifact = 'artifacts_data_raw';
% cfg = ft_rejectartifact(cfg);
% cfg.implicitref   = []; %REF is implicit, what is 0d to, RM is the other one, do average reference with the mastoid

%% ICA
cfg        = [];
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
icomp = ft_componentanalysis(cfg, trials);
%!!! try before vs. after trial segmentation

%% save the ICA unmixing and topolabel, eeg, and eog_bp
icaunmixing = icomp.unmixing;
icatopolabel = icomp.topolabel;

data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '.mat'];
save(data_fname, 'icaunmixing', 'icatopolabel', 'trials', 'eog_trials', 'artifacts_data_raw');


