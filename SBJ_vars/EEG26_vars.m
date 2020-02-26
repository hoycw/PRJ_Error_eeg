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
SBJ_vars.SBJ = 'EEG26';
SBJ_vars.raw_file = {'eeg26.bdf'};
SBJ_vars.bhv_file = 'eeg26_response_log_20200203155436.txt';
SBJ_vars.oddball_file = {'eeg26_oddball_log_20200203153621.txt'};
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
SBJ_vars.ch_lab.bad     = {'PO4', 'P2'};
SBJ_vars.trial_reject_ix = [95, 99, 109, 124, 145, 207, 208, 326, 354, 355, 358, 359, 366, 373, 377, 378, 385, 386, 398, 402, 405, 417, 425, 435, 445, 447, 456, 458, 475, 476, 499, 500, 502, 503, 504, 507, 518, 521, 522, 525, 526, 527, 528, 529, 530, 531, 533, 535, 538, 545, 559, 562, 563, 568, 570, 572, 573, 574, 576, 577, 578, 587, 588];
SBJ_vars.trial_reject_ix_oddball = [11, 110, 133, 179, 180, 207, 212, 236, 261, 360, 373, 374];
SBJ_vars.ica_reject = [2, 5, 10, 13, 16, 19, 20, 24, 25, 27, 32, 33, 34, 35, 36, 40, 44, 47, 51, 52, 56, 57:62];
SBJ_vars.odd_trigger_ix = 3;
SBJ_vars.tt_trigger_ix = 404;
%Messy variance around 800 - 1000


