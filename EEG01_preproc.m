function EEG01_preproc(SBJ, proc_id)
%% Import and preprocess EEG data, then run ICA

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Set up processing and SBJ variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);

%% Load and preprocess the data
% Update ears ref channels with prefix/suffix

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
cfg=[];
cfg.dataset    = SBJ_vars.dirs.raw_filename;
cfg.continuous = 'yes'; 
cfg.lpfilter   = proc_vars.lp_yn;
cfg.hpfilter   = proc_vars.hp_yn;
cfg.bpfilter   = proc_vars.bp_yn;
cfg.bpfreq     = proc_vars.bp_freq;
cfg.bpfiltord = proc_vars.bp_order;
cfg.demean     = proc_vars.demean_yn;
cfg.reref      = proc_vars.reref_yn;
cfg.refmethod  = proc_vars.ref_method;
cfg.refchannel = {ear_lab1, ear_lab2};
data = ft_preprocessing(cfg);

%% Downsample
%Just looked over script a final time! This is one thing I wanted to check
%that its in the right spot that I forgot to mention -- I wanted it to go
%after filtering and trial cutting, but wasn't sure if this was too late!
%(April 13,2019)
if strcmp(proc_vars.resample_yn,'yes')
      cfg=[];
    cfg.resamplefs = proc_vars.resample_freq;
    data = ft_resampledata(cfg, data);
end

%% Fix channel labels
% Remove bad channels
null_neg = {};
for null_ix = 1:numel(SBJ_vars.ch_lab.null)
    null_neg = {null_neg{:},['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.null{null_ix} SBJ_vars.ch_lab.suffix]};
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
cfg.channel = {'all',bad_neg{:},rep_neg{:},ears_neg{:},null_neg{:}};

%cfg.channel = {'all',rep_neg{:},ears_neg{:}};
data = ft_selectdata(cfg,data);

% Strip Pre/Suffix if Necessary
for ch_ix = 1:numel(data.label)
    data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
    data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
    for x = 1:numel(SBJ_vars.ch_lab.replace)
        if strcmp(SBJ_vars.ch_lab.replace{x}{2},data.label{ch_ix})
            data.label{ch_ix} = SBJ_vars.ch_lab.replace{x}{1}; % replaces the label of the externals with the channel they represent
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
% 
% % Add in trigger channel to segment trials
% cfg.channel = SBJ_vars.ch_lab.trigger;
% trigger = ft_selectdata(cfg,data);

% Combine bipolar EOG
cfg = [];
eog = ft_appenddata(cfg,eog_h,eog_v);%,trigger

% Remove unipolar EOG
warning('WARNING!!! Assuming Fp2 is the second vertical EOG!');
eog_v_low_ix = ~strcmp(SBJ_vars.ch_lab.eog_v,'Fp2');    % only toss the lower one
eog_neg = [fn_ch_lab_negate(SBJ_vars.ch_lab.eog_h),fn_ch_lab_negate(SBJ_vars.ch_lab.eog_v(eog_v_low_ix))];
trig_neg = fn_ch_lab_negate({SBJ_vars.ch_lab.trigger});
cfg = [];
cfg.channel = {'all',eog_neg{:},trig_neg{1}};
data = ft_selectdata(cfg,data);

% Baseline-correction options- Not doing right now
%cfg.demean          = 'yes'; %baseline correct data before you resample?
%cfg.baselinewindow  = [-0.25 -0.05]; %set baseline to -.2 to 0 seconds, then correct it before you resample


%!!! WTF was this about?
% cfg = [];
% cfg.latency = [-0.2 1];
% trials = ft_selectdata(cfg,trials);

% cfg.artfctdef.visual.artifact = 'artifacts_data_raw';
% cfg = ft_rejectartifact(cfg);
% cfg.implicitref   = []; %REF is implicit, what is 0d to, RM is the other one, do average reference with the mastoid

% %% Cut into trials
% % Must segment before downsampling because trigger channel read from
% % original file
% % Need to cut into trials here so that can NaN out bad trials later in
% % script
% cfg = [];
% cfg.dataset             = SBJ_vars.dirs.raw_filename;
% cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
% cfg.trialdef.eventvalue = proc_vars.event_code;        % feedback cocde
% cfg.trialdef.prestim    = proc_vars.trial_lim_s(1);
% cfg.trialdef.poststim   = proc_vars.trial_lim_s(2);
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
load([SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat']);
cfgs        = [];
cfgs.hpfilter = proc_vars.ICA_hp_yn
cfgs.hpfreq = proc_vars.ICA_hp_freq
cfgs.artfctdef.visual.artifact = bad_epochs;
cfgs.artfctdef.reject = 'nan';
data = ft_rejectartifact(cfgs, data);
%cfgs.trials.trial(bad_raw_trials) = NaN;
data_ICA   = ft_preprocessing(cfgs, data)
cfgs.method = 'runica'; % this is the default and uses the implementation from EEGLAB
%cfgs.channel = setxor(data.label(1:64), SBJ_vars.ch_lab.bad);
icomp = ft_componentanalysis(cfg, data_ICA);

%!!! try before vs. after trial segmentation

%% save the ICA unmixing and topolabel, eeg, and eog_bp
icaunmixing = icomp.unmixing;
icatopolabel = icomp.topolabel;

data_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
save(data_fname, 'icaunmixing', 'icatopolabel', 'data', 'eog','bad_epochs');% 'trials', 'eog_trials',


