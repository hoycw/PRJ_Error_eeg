%% EEG Pilot 06 Processing Variables
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EP06';
SBJ_vars.raw_file = 'TT_Cyclone_pilot06.bdf';
SBJ_vars.bhv_file = 'TT_Cyclone_pilot06_response_log_20180426140027.txt';
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
SBJ_vars.ch_lab.replace = {{'F6', 'EXG6'},{'P8','EXG7'},{'P2','EXG8'}}; % {{'final','EXG#'},{'final2','EXG#2'}}
SBJ_vars.ch_lab.prefix  = '1-';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {...
    'Iz','Oz','PO8','POz','PO3','O2',... % noisy channels
    'T8','F6','AF3','TP7','PO7','O1'... % noisy channels
    };
%SBJ_vars.ref_exclude = {}; %exclude from the CAR

%SBJ_vars.trial_reject_ix = [...
    %142, 228, 375, 376, 387, 398, 402, 409,...
    %410, 420, 441, 462, 463, 479, 490, 497,...
    %513, 533, 543, 548, 552:554, 570, 580];
SBJ_vars.trial_reject_n = [141, 227, 374, 375, 386, 397, 401, 408, 409, 419, 440, 461, 462, 478, 489, 496, 512, 532, 542, 547, 551, 552, 553, 569, 579];
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
    % POz has big ripples, PO3 at times too
    % TP8 also has big fluctuations
    % PO8 has big noise at times
    % O2 breaks loose at some point
    % AF3 becomes very onisy at some point
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

%--------------------------------------
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
