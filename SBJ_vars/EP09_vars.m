%% EEG Pilot 09 Processing Variables
% Battery died, EEG restarted but not paradigm
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EP09';
SBJ_vars.raw_file = {'pilot09-2.bdf'};
SBJ_vars.bhv_file = 'pilot09_response_log_20180426140027.txt';
SBJ_vars.block_name = {'r1'};
%NOTE: I ended up ignoring the whole first set because in my notes it says
%that the participent was misunderstanding the task for all of block 1 and
%hte battery cut out at the beginning of block 2.
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
SBJ_vars.ch_lab.replace = {}; % {{'final','EXG#'},{'final2','EXG#2'}}
SBJ_vars.ch_lab.prefix  = '1-';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {'P2', 'P9', 'P1'};
SBJ_vars.ch_lab.null    = {'EXG6', 'EXG7', 'EXG8'};
SBJ_vars.tt_trigger_ix = 3;
SBJ_vars.odd_trigger_ix = 10000;
SBJ_vars.trial_reject_ix = [1:2, 6:8, 11, 17, 20, 24, 27, 35, 42, 45, 47, 50, 64, 74:76, 104, 135, 137, 139, 140, 171, 201, 238:243, 273, 282, 283, 286, 310, 311, 312, 313, 316, 317, 323, 324, 325, 335, 337, 349, 402:406, 408:409, 416];
SBJ_vars.ica_reject = [1 2 6 8 9 17 15 16 17 12 21 25 36 30 37 39 44 43 59 51 53 47 48 58 60];
