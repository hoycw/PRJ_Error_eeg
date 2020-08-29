%% EEG 04 Processing Variables
% EEG24 restarted paradigm and data recording (battery died)
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EEG24';
SBJ_vars.raw_file = {'eeg24.bdf', 'eeg24_2.bdf'};
SBJ_vars.bhv_file = {'eeg24_response_log_20191208101056.txt', 'eeg24_2_response_log_20191208103858.txt'};
SBJ_vars.oddball_file = 'eeg24_oddball_log_20191208095654.txt';
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
SBJ_vars.ch_lab.bad     = {'T8', 'FT7', 'AF8', 'F6'
    };
SBJ_vars.trial_reject_ix = [56, 75, 76, 92, 126, 148, 221, 222, 271, 294, 323, 324, 339, 340, 157, 360, 428, 429, 468, 562];
SBJ_vars.trial_reject_ix_oddball = [1, 30, 141, 149, 164, 169, 171, 224, 260, 261, 279, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 295, 296, 297, 298, 300, 301, 302, 303, 316, 344, 346, 359, 360, 363, 364, 367, 370, 372, 373, 374, 375, 378, 379];
SBJ_vars.ica_reject = [1, 5, 7, 9, 10, 11, 12, 16, 17, 23, 25, 26, 28, 29, 31, 34, 36:40, 43:60];
SBJ_vars.tt_trigger_ix = 3;
SBJ_vars.odd_trigger_ix = 404;
SBJ_vars.endsample = {1114624};
SBJ_vars.origsample = 512;
SBJ_vars.resample = 250;
