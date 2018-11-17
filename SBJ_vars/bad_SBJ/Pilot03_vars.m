%% IR39 Processing Variables
addpath('/Users/SCS22/Documents/MATLAB/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'Pilot03';
SBJ_vars.raw_file = 'Pilot03.bdf' ;
SBJ_vars.block_prefix = '';

SBJ_vars.dirs.SBJ     = ['/Users/SCS22/Desktop/Preprocessing_Work/' SBJ_vars.SBJ '/'];
%SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
%SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
%SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
%SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
%SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
%if ~exist(SBJ_vars.dirs.import,'dir')
    %mkdir(SBJ_vars.dirs.import);
%end
%if ~exist(SBJ_vars.dirs.preproc,'dir')
 %   mkdir(SBJ_vars.dirs.preproc);
%end
%if ~exist(SBJ_vars.dirs.events,'dir')
%    mkdir(SBJ_vars.dirs.events);
%end
%if ~exist(SBJ_vars.dirs.proc,'dir')
%    mkdir(SBJ_vars.dirs.proc);
%end

SBJ_vars.dirs.raw_filename = strcat('/Users/SCS22/Desktop/Preprocessing_Work/', SBJ_vars.SBJ, '/', SBJ_vars.raw_file);

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

SBJ_vars.ears = {'1-EXG1', '1-EXG2'};
SBJ_vars.horiz = {'EXG3', 'EXG4'};
SBJ_vars.vert = {'EXG5', 'Fp2'};
SBJ_vars.replacement_channels = {{}, {}};

%SBJ_vars.ch_lab.bad = {...
    %'LAM*','LHH1','LHH2','LHH3','LTH1','LTH2','LTH3','LTH10',...%Epileptic
    %'RHH1','RHH2','RHH3','RTH2','RTH3','RTH4','RAM1','RAM2','RAM3','RAM6',...% Epileptic
    %'RIN10','ROF10','LIN3','LOF1','LOF10','LIN10','ROF10','LHH10',...%noisy or out of brain
    %'PMF*','AMF*','OF*','AT*','PT*','AD*','HD*','DC03','DC04',....% Not real data
    %'LT1','LT2','LT3','LT4','LT5','LT6',...% Not real data, using 'LT*' tosses LTH elecs too
    %'E','V1',...% Not real data
    %'EKG*'...
    %};
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
SBJ_vars.analysis_time = {[660000:1662036]};%660000/512 sample #/sampling frequency??? still not 100% positive on this

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
SBJ_vars.trial_reject_ix = [];
SBJ_vars.channels_reject_ix = {};
%--------------------------------------
% Component Paramaters
%--------------------------------------
SBJ_vars.top_comp_cut = 0.1;
SBJ_vars.rejcomp = [];
