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
SBJ_vars.SBJ = 'EEG17';
SBJ_vars.raw_file = {'eeg17.bdf'};
SBJ_vars.bhv_file = {'eeg17_response_log_20191106172905.txt'};
SBJ_vars.oddball_file = 'eeg17_oddball_log_20191106171454.txt';
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
SBJ_vars.ch_lab.null = {'EXG6', 'EXG7'};
SBJ_vars.ch_lab.replace = {}; % {{'final','EXG#'},{'final2','EXG#2'}}
SBJ_vars.ch_lab.prefix  = '';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {'T7', 'Iz'
    };
%Maybe get rid of AF7
SBJ_vars.trial_reject_ix = [22, 33, 40, 99, 144, 190, 233, 234, 275, 276, 284, 285, 315, 316, 335, 336, 470, 485, 491, 518, 521, 524, 582];
SBJ_vars.trial_reject_ix_oddball = [121, 207, 208,  324, 325];
SBJ_vars.ica_reject = [1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 14, 16, 17, 19, 20, 21, 22, 28, 29, 34, 32, 35, 36, 38, 41, 44, 45, 49, 50, 51, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63];
% Maybe too many components!!
SBJ_vars.tt_trigger_ix = 3;
SBJ_vars.odd_trigger_ix = 404;
%Varience = 200 - 600, higher end
