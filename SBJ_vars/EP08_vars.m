%% Pilot08 Processing Variables
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EP08';
SBJ_vars.raw_file = 'Pilot08.bdf';
SBJ_vars.bhv_file = 'Pilot08_response_log_20181101084314.txt';
SBJ_vars.block_prefix = '';

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Error_eeg/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
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
SBJ_vars.ch_lab.replace = {{'O2', 'EXG6'}}; % {{'final','EXG#'},{'final2','EXG#2'}}
SBJ_vars.ch_lab.prefix  = '1-';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {'T8', 'FT8', 'PO4', 'C6', 'PO8', 'F8', 'F6', 'AF8', 'O2', 'P2','Iz','Pz'};
%SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.trial_reject_n = [503, 504, 406];
SBJ_vars.ica_reject = [1 7 10 13 18 19 23 26 27 28 25 30 44 40 46 50 55];
%SBJ_vars.trial_reject_ix = [88, 94, 103, 154, 157, 217, 327, 343, 370, 406, 417, 440, 490, 499, 511, 545];
%SBJ_vars.trial_reject_n = [87, 93, 102, 153, 156, 216, 326, 342, 369, 405, 416, 439, 489, 510, 544];
%--------------------------------------
% Noise Notes
%--------------------------------------
% recording info sheet notes:
    %'PO7','O1','Iz','Oz'... % noisy channels
% PSD Notes:
    %'AF3','CP2','Iz','Oz','O1','O2','PO3','PO7','PO8','POz' - PSD looks noisy
    %'F1' - strange flat PSD
    %'F6','P2','P8' empty
% Raw View notes:
    % FT8 spiking, toss it
    % T8 messy, f8 and f6 messy, af8 messy
    % PO8 has big noise at times
    % FT8 loose at 170- 200 seconds, will get rid of anyways
    % P2 gets weird at 230s
    %IZ gets messy 970
    %02 1680-1700
% databrowser post-IC rejection:
    % channels: Iz, Oz, PO8, PO3, T8 (trial 31), POz (esp. t 151), 216 starts O2 loose, F6 (t 351), AF3 (t 413), F4(t 417), TP7 (t 436)
    % trials: 142 (EOG missed?), 228, 375, 376, 387, 398 (P5), 402, 409, 410, 420, 441, 462, 463, 479, 490, 497, 513, 533, 543, 548, 552:554, 559?, 570, 576?, 580
% ft summary:
    % Fp1, AF7, AF3, PO3, Iz, Oz, POz, Fpz, Fp2, AF8, F4, T8, PO8, O2, F6
    % trials (n/582): 31, 32, 329

% pre-ICA rejection:
% ft summary notes:
    % channels (n/64): Iz, Oz, PO8, O2
    % trials (n/582): 151, 351, 386, 413, 419, 462, 463, 513, 554
% ft summary EOG notes:
    % trials (n/582): 139, 142, 325, 542, 554, 559, 576, 580bp

%--------------------------------------
% Time Parameters
%--------------------------------------
% SBJ_vars.analysis_time = {};

%-------------------------s-------------
% Trials and Channels to Reject
%--------------------------------------
% These should be indices AFTER SBJ05 has run!
% original trial_reject_ix = [52 61 77 125 154 156 185 187 131 132 133 205 254 265 280 283 4303 311 318 319 320];
%SBJ_vars.trial_reject_ix = [52 61 77 125 154 156 185 187 131 132 133 205 254 265 280 283 303 311 318 319 320];
%SBJ_vars.channels_reject_ix = {'T7','T8'};

%--------------------------------------
% Component Paramaters
%--------------------------------------
% SBJ_vars.top_comp_cut = 0.1;
%SBJ_vars.rejcomp = [1];
