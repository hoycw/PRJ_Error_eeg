function clean_data = fn_load_clean_experiment(SBJ, proc_id)
%% Load, clean, and reconstruct full experiment time series
%   Essentially reproduces SBJ02b but without trial segmentation
% INPUTS:
%   SBJ [str] - ID of subject
%   proc_id [str] - ID of preprocessing pipeline
% OUTPUTS:
%   clean_data [FT struct] - data after rejecting ICs and repairing channels

%% Load Subject Data
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Reload Data
load([SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'],'data','icaunmixing','icatopolabel');
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat'], 'heog_ics','veog_ics');

%% Rebuild Independent Components from Full Session Data
cfg           = [];
cfg.unmixing  = icaunmixing;
cfg.topolabel = icatopolabel;
ica           = ft_componentanalysis(cfg, data);

%% Rebuild Data with Only Clean ICA Components
cfg = [];
cfg.component = unique([SBJ_vars.ica_reject, heog_ics, veog_ics]);
cfg.demean = 'no';
clean_data = ft_rejectcomponent(cfg, ica);

%% Repair Bad Channels
cfg = [];
cfg.method         = 'average';
cfg.missingchannel = SBJ_vars.ch_lab.bad(:); % not in data (excluded from ica)
cfg.layout         = 'biosemi64.lay';

% Identify spatial relationships between neighboring channels
cfgn = [];
cfgn.channel = 'all';
cfgn.layout  = 'biosemi64.lay';
cfgn.method  = 'template';
cfg.neighbours = ft_prepare_neighbours(cfgn);

clean_data = ft_channelrepair(cfg, clean_data);

end

