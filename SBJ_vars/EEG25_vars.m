%% EEG 00 Processing Variables
% EEG25 restarted paradigm and data recording (battery died)
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EEG25';
SBJ_vars.raw_file = {'eeg25.bdf', 'eeg25block8.bdf'};
SBJ_vars.bhv_file = {'eeg25_response_log_20200126163741.txt', 'eeg25_block 8_response_log_20200126172226.txt'};
SBJ_vars.oddball_file = {'eeg25_oddball_log_20200126162425.txt'};
SBJ_vars.block_name = {'1', '2'};
SBJ_vars.odd_block_name = {''};

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Error_eeg/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
SBJ_vars.dirs.proc_stack    = [SBJ_vars.dirs.SBJ '04_proc/plot/'];
if ~exist(SBJ_vars.dirs.import,'dir')
mkdir(SBJ_vars.dirs.import);
end
if ~exist(SBJ_vars.dirs.preproc,'dir')
mkdir(SBJ_vars.dirs.preproc);
end
if ~exist(SBJ_vars.dirs.events,'dir')
mkdir(SBJ_vars.dirs.events);
end
if ~exist(SBJ_vars.dirs.proc,'dir')
mkdir(SBJ_vars.dirs.proc);
end

SBJ_vars.dirs.raw_filename = strcat(SBJ_vars.dirs.raw, SBJ_vars.raw_file);

%--------------------------------------
% Channel Selection
%--------------------------------------
% Channel Labels
SBJ_vars.ch_lab.ears    = {'EXG1', 'EXG2'};
SBJ_vars.ch_lab.eog_h   = {'EXG3', 'EXG4'};
SBJ_vars.ch_lab.eog_v   = {'EXG5', 'Fp2'};
SBJ_vars.ch_lab.null = {'EXG6', 'EXG7', 'EXG8'};
SBJ_vars.ch_lab.replace = {}; % {{'final','EXG#'},{'final2','EXG#2'}}
SBJ_vars.ch_lab.prefix  = '';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {'T7', 'T8'};
SBJ_vars.trial_reject_ix = [18, 317];
SBJ_vars.trial_reject_ix_oddball = [29, 103, 242, 374];
SBJ_vars.ica_reject = [2, 7, 8, 9, 14, 19, 20, 21, 24, 25, 26, 38, 41, 45, 51, 54:56, 59, 60, 61];
SBJ_vars.odd_trigger_ix = 404;
SBJ_vars.tt_trigger_ix = 3;
SBJ_vars.endsample = {259072}; %1620992
SBJ_vars.resample = 250;
SBJ_vars.origsample = 512;
