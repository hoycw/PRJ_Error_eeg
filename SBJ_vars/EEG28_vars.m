%% EEG 00 Processing Variables
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('/Users/aasthashah/','dir'); root_dir='/Users/aasthashah/Desktop/';ft_dir='/Users/aasthashah/Desktop/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EEG28';
SBJ_vars.raw_file = {'eeg28.bdf'};
SBJ_vars.bhv_file = 'eeg28_response_log_20200205154520.txt';
SBJ_vars.oddball_file = {'eeg28_oddball_log_20200205151020.txt'};
SBJ_vars.block_name = {''};

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
SBJ_vars.ch_lab.bad     = {'PO7', 'P2', 'P8', 'CP3', 'TP8', 'P10'};
SBJ_vars.trial_reject_ix = [55, 74, 78, 85, 87, 140, 141, 154, 174, 191, 223, 251, 260, 271, 286, 329, 351, 365, 369, 377, 382, 394, 404, 433, 442, 478, 522, 533, 536, 546, 548, 579];
SBJ_vars.trial_reject_ix_oddball = [12, 80, 134];
SBJ_vars.ica_reject = [1, 4, 5, 6, 8, 9, 11, 15, 16, 18, 30, 34, 36, 38, 39, 40, 41, 48, 49, 50, 52:58]; %26??
SBJ_vars.odd_trigger_ix = 3;
SBJ_vars.tt_trigger_ix = 404;
%This data set may be too messy! Crazy blinks
%Max variance around 800 now
%Max diff around 8


