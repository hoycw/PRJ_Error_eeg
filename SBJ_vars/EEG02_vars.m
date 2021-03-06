%% EEG 02 Processing Variables
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EEG02';
SBJ_vars.raw_file = {'eeg02.bdf'};
SBJ_vars.bhv_file = 'eeg02_response_log_20190715153326.txt';
SBJ_vars.oddball_file = 'eeg02oddball_oddball_log_20190715151729.txt';
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
SBJ_vars.ch_lab.prefix  = '1-';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {'Iz'};
SBJ_vars.trial_reject_ix = [11, 12, 13, 34, 68, 73, 81, 93, 99, 100, 69, 101, 103, 104, 107, 109, 112, 125, 117, 120, 123, 124, 131, 159, 160, 189, 197, 205, 206, 225, 226, 238, 256, 269, 273, 282, 283, 287, 297, 307, 314, 316, 335, 336, 338, 339, 359, 372, 373, 374, 384, 311, 399, 439, 445, 532];
SBJ_vars.ica_reject = [1, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 17, 21, 29, 35, 38, 39, 54, 56, 57, 58,63];
%Maybe 6
SBJ_vars.trial_reject_ix_oddball = [50, 160, 370, 149, 25, 37, 43, 50, 51, 86, 99, 102, 149, 160, 312, 309, 326, 381, 344]; 
SBJ_vars.tt_trigger_ix = 3;
SBJ_vars.odd_trigger_ix = 403;
% Variance is really high, lots of bad trials
% Max: Around 30
%REJECT
