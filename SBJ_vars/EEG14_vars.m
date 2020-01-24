%% EEG 04 Processing Variables
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
SBJ_vars.SBJ = 'EEG14';
SBJ_vars.raw_file = {'eeg14.bdf'};
SBJ_vars.bhv_file = {'eeg14_response_log_20191008170025.txt'};
SBJ_vars.oddball_file = 'eeg14_oddball_log_20191008164516.txt';
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
SBJ_vars.ch_lab.bad     = {'POz', 'P2'}
SBJ_vars.trial_reject_ix = [129, 154, 181, 241, 242, 297, 305, 316, 317, 318, 319, 331, 338, 339, 340, 363, 364, 375, 379, 428, 429, 430, 448, 467, 469];
SBJ_vars.trial_reject_ix_oddball = [133, 169, 198, 292, 323, 330, 359, 361, 362, 380, 381];
SBJ_vars.ica_reject = [1, 4, 9, 11, 11, 14, 16, 19, 20, 21, 23, 24, 26, 28, 30, 31, 32, 33, 36, 37, 42, 43, 46, 48, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62];
SBJ_vars.tt_trigger_ix = 3;
SBJ_vars.odd_trigger_ix = 404;
% Max Variance: 800, most below 200 (after rejected) -- pretty messy got
% rid of a lot of trials
% Max is around 40
%Messy data set!!!
%Variance is 100 - 600
