%% EEG 12 Processing Variables
% EEG12 restarted paradigm
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EEG12';
SBJ_vars.raw_file = {'eeg12.bdf'};
SBJ_vars.bhv_file = {'eeg12(part 1)_response_log_2019081914'};
SBJ_vars.oddball_file = {'eeg12_oddball_log_20190819141030.txt','eeg12_oddball_log_20190819141120.txt','eeg12_oddball_log_20190819142106.txt'};
SBJ_vars.block_name = {''};
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
SBJ_vars.ch_lab.null = {'EXG8'};
SBJ_vars.ch_lab.replace = {};
% {{'P2','EXG6'},{'PO4','EXG7'}}; % {{'final','EXG#'},{'final2','EXG#2'}}
% THESE ARE BAD
SBJ_vars.ch_lab.prefix  = '1-';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {...
'T7','P9','AF8', 'EXG6', 'EXG7'
    };
SBJ_vars.trial_reject_ix = [53 54 55 56 57 73 74 105 106 118 123 140 144 167 174 175 189 190 195 209 210 211 212 213 216 244 260 267 268 278 283 309 329 330 388 389 411 412 433 460 467 493 496 497 503 504 505 525 535 568];
%old: SBJ_vars.trial_reject_ix = [21 42 71 100 101 128 129 183 209 212 213 214 215 219 149 220 221 242 252];
%old: SBJ_vars.trial_reject_ix_oddball = [25 91 115 116 121 158 276 279 280 281 292 293 297 308 309 310 311 344 345 346 347 363 375 376 377 378];
SBJ_vars.trial_reject_ix_oddball = [113 114 115 116 121 279 280 292 293 297 309 311 344 345 346 347 363 375 376 377 390];
SBJ_vars.ica_reject = [1 3 4 6 8 9 10 11 13 16 23 25 26 28 38 43 44 48 51 53 54:61];
SBJ_vars.tt_trigger_ix = 472;
SBJ_vars.odd_trigger_ix = 3;
% Variance highish -- around 700
%Max around 14

