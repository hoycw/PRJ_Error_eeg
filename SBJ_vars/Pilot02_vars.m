%% Pilot02 Processing Variables
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath(ft_dir);
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'Pilot02';
SBJ_vars.raw_file = 'Pilot02.bdf' ;
SBJ_vars.block_prefix = '';

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Error_eeg/data/' SBJ_vars.SBJ '/'];
data_out_filename_ICA = [SBJ_vars.dirs.SBJ 'Pilot02ICASorted' datestr(now,'mm-dd-yyyy HH-MM')];
data_out_filename_cleaned = [SBJ_vars.dirs.SBJ 'Pilot02clean' datestr(now,'mm-dd-yyyy HH-MM')];
data_out_filename_erp = [SBJ_vars.dirs.SBJ  'erps' datestr(now,'mm-dd-yyyy HH-MM')]
data_out_filename_TFA = [SBJ_vars.dirs.SBJ 'TFA' datestr(now,'mm-dd-yyyy HH-MM')]
data_out_filename_rawdatainspect = [SBJ_vars.dirs.SBJ 'rawdatainspect' datestr(now,'mm-dd-yyyy HH-MM')]
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
hdr = ft_read_header(SBJ_vars.dirs.raw_filename);
SBJ_vars.orig_n_ch = length(hdr.label);
SBJ_vars.orig_n_samples = hdr.nSamples;
SBJ_vars.orig_srate = hdr.Fs;
clear hdr;


%SBJ_vars.ch_lab.probes = {'ROF','RAC','RIN','RAM','RHH','RTH','LOF','LIN','LHH','LTH'};%'LAM' is all bad
%SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};%'BP'
%SBJ_vars.ch_lab.ROI    = {'RAC*','ROF*','RIN*','LOF*','LIN*'};

SBJ_vars.ch_lab.ears = {'EXG1', 'EXG2'};
SBJ_vars.ch_lab.eog_h = {'EXG3', 'EXG4'};
SBJ_vars.ch_lab.eog_v = {'EXG5', 'Fp2'};
SBJ_vars.ch_lab.replacements = {{'Fp1', 'EXG6'}}; % {{'final','EXG#'},{'final2','EXG#2'}}
SBJ_vars.ch_lab.prefix = '1-';    % before every channel
SBJ_vars.ch_lab.suffix = '';    % after every channel

SBJ_vars.ch_lab.bad = {...
    'EXG7','EXG8'... % empty external channels
    };
%SBJ_vars.ref_exclude = {}; %exclude from the CAR
%SBJ_vars.ch_lab.eeg = {'C3','CZ','C4','FZ','OZ'};
%SBJ_vars.ch_lab.eeg_ROI = {'CZ','FZ'};
%SBJ_vars.ch_lab.lap_ref = {{'C3','C4'},{'C3','C4'}};
%SBJ_vars.ch_lab.eog = {'LSH','LLE','RSH'}; % lower left, upper right, ???
%SBJ_vars.ch_lab.photod = {'DC01'};
%SBJ_vars.ch_lab.mic    = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% most have only regular harmonics with normal width, and 120 and 240 are weak
% RBT, RPIN, RSMA,LAC have an extra peak at 200
% RHH6 has really bad at all harmonics, like LUE and other nonsense chan
%SBJ_vars.notch_freqs = [60 120 180 240 300]; %200 shoudl come out in re-referencing
%SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% events start ~155 or 160s to ~384; ~580 to end (~1360?)
SBJ_vars.analysis_time = {};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
%SBJ_vars.artifact_params.std_limit_raw = 7;
%SBJ_vars.artifact_params.hard_threshold_raw = 400;

%SBJ_vars.artifact_params.std_limit_diff = 7;
%SBJ_vars.artifact_params.hard_threshold_diff = 40;

%--------------------------------------
% Trials and Channels to Reject
%--------------------------------------
% These should be indices AFTER SBJ05 has run!
SBJ_vars.trial_reject_ix = [52 61 77 125 154 156 185 187 131 132 133 205 254 265 280 283 4303 311 318 319 320];
SBJ_vars.channels_reject_ix = {'T7','T8'};
%--------------------------------------
% Component Paramaters
%--------------------------------------
SBJ_vars.top_comp_cut = 0.1;
SBJ_vars.rejcomp = [1];
