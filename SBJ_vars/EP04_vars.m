%% Pilot04 Processing Variables
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'EP04';
SBJ_vars.raw_file = 'Pilot04.bdf' ;
SBJ_vars.bhv_file = 'Pilot4_response_log_20180322141434.txt' ;
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
SBJ_vars.ch_lab.replace = {{'Fp1', 'EXG6'}}; % {{'final','EXG#'},{'final2','EXG#2'}}
SBJ_vars.ch_lab.prefix  = '1-';    % before every channel
SBJ_vars.ch_lab.suffix  = '';    % after every channel
SBJ_vars.ch_lab.trigger = 'Status';
SBJ_vars.ch_lab.bad     = {...
    'EXG7','EXG8'... % empty external channels
    };
%SBJ_vars.ref_exclude = {}; %exclude from the CAR

%--------------------------------------
% Time Parameters
%--------------------------------------
% events start ~155 or 160s to ~384; ~580 to end (~1360?)
SBJ_vars.analysis_time = {};

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
