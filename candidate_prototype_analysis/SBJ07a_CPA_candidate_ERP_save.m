function SBJ07a_CPA_candidate_ERP_save(SBJ,eeg_proc_id,odd_proc_id,an_id,cpa_id)
%% Plot ERPs for single SBJ
% INPUTS:
%   conditions [str] - group of condition labels to segregate trials

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' eeg_proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cpa_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' cpa_id '_vars.m'];
eval(cpa_vars_cmd);

% Load data
load([SBJ_vars.dirs.preproc SBJ '_' eeg_proc_id '_02a.mat'],'ica');
load([SBJ_vars.dirs.proc SBJ '_' cpa_id '_' odd_proc_id '_prototype.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' eeg_proc_id '_final.mat']);

%% Reconstruct and Clean Data
% Reconstruction
cfg = [];
cfg.component = setdiff(1:numel(ica.label),final_ics);
recon = ft_rejectcomponent(cfg, ica);

% Repair Bad Channels
cfg = [];
cfg.method         = 'average';
cfg.missingchannel = SBJ_vars.ch_lab.bad(:); % not in data (excluded from ica)
cfg.layout         = 'biosemi64.lay';

cfgn = [];
cfgn.channel = 'all';
cfgn.layout  = 'biosemi64.lay';
cfgn.method  = 'template';
cfg.neighbours = ft_prepare_neighbours(cfgn);

full_recon = ft_channelrepair(cfg, recon);

% Check for null channels (EEG17-23, EEG29)
null_neg = {};
for null_ix = 1:numel(SBJ_vars.ch_lab.null)
    null_lab = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.null{null_ix} SBJ_vars.ch_lab.suffix];
    if any(strcmp(full_recon.label,null_lab))
        null_neg = [null_neg {['-' null_lab]}];
    end
end

% Toss bad trials
cfgs = [];
cfgs.channel = [{'all'}, null_neg];
cfgs.trials  = setdiff([1:numel(full_recon.trial)], SBJ_vars.trial_reject_ix);
clean_trials = ft_selectdata(cfgs, full_recon);

%% Select trials and epoch for plotting
% Select epoch
% Realign data to desired event
if ~strcmp(proc.event_type,an.event_type)
    cfg = [];
    % Match desired time to closest sample index
    if strcmp(proc.event_type,'S') && strcmp(an.event_type,'F')
        prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);
        cfg.offset = -(prdm_vars.target + prdm_vars.fb_delay)*clean_trials.fsample;
    elseif strcmp(proc.event_type,'S') && strcmp(an.event_type,'R')
        cfg.offset = round(-bhv.rt*clean_trials.fsample);
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

% Check window consistency
%   Check trial_lim_s is within trial time (round to avoid annoying computer math)
if round(an.trial_lim_s(1)+1/roi.fsample,3) < round(roi.time{1}(1),3) || ...
        round(an.trial_lim_s(2)-1/roi.fsample,3) > round(roi.time{1}(end),3)
    error('an.trial_lim_s is outside data time bounds!');
end

% Select window and channels of interest
cfgs = [];
cfgs.channel = an.ROI;
cfgs.latency = an.trial_lim_s;
roi = ft_selectdata(cfgs, roi);

%% Preprocess Data for ERP
cfgpp = [];
cfgpp.hpfilter       = an.hp_yn;
cfgpp.hpfreq         = an.hp_freq;
cfgpp.hpfiltord      = an.hp_filtord; % Leaving blank causes instability error, 1 or 2 works 
cfgpp.lpfilter       = an.lp_yn;
cfgpp.lpfreq         = an.lp_freq;
cfgpp.demean         = an.demean_yn;
cfgpp.baselinewindow = an.bsln_lim;
if isfield(an,'hilbert')
    cfgpp.hilbert = an.hilbert;
end
roi = ft_preprocessing(cfgpp, roi);

%% Downsample
if an.dsamp_yn
    cfgds = [];
    cfgds.resamplefs = an.dsamp_freq;
    cfgds.detrend    = 'no';
    roi = ft_resampledata(cfgds, roi);
end

%% Save Results
data_out_fname = [SBJ_vars.dirs.proc SBJ '_' cpa_id '_' an_id '.mat'];
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','roi');

end
