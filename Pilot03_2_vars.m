%% preprocessing --ERP practice
%data_dir = '/Users/SCS22/Desktop/Pilot2_ICATesting/pilot01/'
%data_filename = [data_dir 'Pilot01.bdf'];
%data_out_filename = [data_dir 'ICAsorted_pilot1_take2'];
%top_comp_cut = 0.1;
%replacement_channels = {{'1-CP2', '1-P2', '1-P8'}, {'1-EX2', '1-EX3', '1-EX6'}};
%ears = {'1-EXG1', '1-EXG2'};
%horiz = {'1-EXG3', '1-EXG4'};
%vert = {'1-EXG5', '1-Fp2'};
%addpath('/Users/SCS22/Documents/MATLAB/fieldtrip/');
%ft_defaults
function SBJ01_import_data(SBJ)

addpath(genpath('/Users/SCS22/Desktop/Knight_Lab/Preprocessing_Work/'));
addpath('/Users/SCS22/Documents/MATLAB/fieldtrip/');
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
load(data_in_filename, 'browsed_data_raw', 'artifacts_data_raw','rawdata')
%% referencing
  % L and R ear lobe (didn't record mastoids)

% Baseline-correction options- Not doing right now
%cfg.demean          = 'yes'; %baseline correct data before you resample?
%cfg.baselinewindow  = [-0.25 -0.05]; %set baseline to -.2 to 0 seconds, then correct it before you resample
cfg.artfctdef.visual.artifact = artifacts_dataraw;
cleandata = ft_rejectartifact(cfg, rawdata);

%% Prepare to Import Data
% eval(strcat(scripts_dir,'SBJName_vars.m'));
cfg = [];
cfg.dataset             = cleandata
cfg.reref = 'yes';
cfg.implicitref   = []; %REF is implicit, what is 0d to, RM is the other one, do average reference with the mastoid
cfg.refchannel    = SBJ_vars.ears; 
% Filtering
cfg.lpfilter = 'no';
cfg.hpfilter = 'no';
cfg.bpfilter = 'yes';
%low pass filter at 8 hz
cfg.bpfreq   = [0.1 200]; 
% Define trials to prepare to cut into trials
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = 2;
cfg.trialdef.prestim    = -.3;
cfg.trialdef.poststim   = 1;
cfg.trialfun            = 'tt_trialfun';
cfg = ft_definetrial(cfg);

data = ft_preprocessing(cfg);


%% Resampling
cfg=[];
cfg.resamplefs = 250;
data = ft_resampledata(cfg, data);

cfg = [];
badindex = match_str(data.label, SBJ_vars.replacement_channels{1});
cfg.channel = setdiff([1:numel(data.label)], badindex);
%cfg.channel = setdiff(data.label,replacementchannels{1});
data = ft_selectdata(cfg, data);
%try to remove the replacement channels here by changing channel, wasnt
%workign with strings
%% Get rid of 1s prefix
for i = 1: numel(data.label)
   data.label{i} = data.label{i}(3:end);
   for x = 1: numel(SBJ_vars.replacement_channels{1})
        if (strcmp(SBJ_vars.replacement_channels{2}{x},data.label{i}))
           data.label{i} = SBJ_vars.replacement_channels{1}{x}; % replaces the label of the externals with the channel they represent
        end
   end
end
%% Extract EOGs and Bipolar Rereference (based off Arjen's)
% Extract EOG channels
cfg = [];
cfg.channel = ft_channelselection({'EXG*', 'Fp2'}, data.label);  % before was using ft_channelselection
eog = ft_selectdata(cfg, data);

%% Prepare the bipolar version
cfg.channel = SBJ_vars.ears; %%had some trouble with this line working too
eog_bp = eog;
eog_bp.label{1} = 'horz';
eog_bp.label{2} = 'vert';

for t = 1:length(eog.trial)
    %!!! get rid fo the ear channels!
    eog_bp.trial{1, t}(1, :) = eog.trial{1,t}(5, :) - eog.trial{1,t}(4, :); %% again, how to select strings instead of indexing into cell array?? to pick out horiz(1) - horiz(2)
    eog_bp.trial{1, t}(2, :) = eog.trial{1,t}(6, :) - eog.trial{1,t}(1, :); %vertical
end

%% Select EEG Data
cfg = [];
cfg.channel = {'all', '-EXG*', '-atus'}; % leave in FP2
eeg = ft_selectdata(cfg, data);

%% ICA
cfg        = [];
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
icomp = ft_componentanalysis(cfg, eeg);

%% save the ICA unmixing and topolabel, eeg, and eog_bp
icaunmixing = icomp.unmixing;
icatopolabel = icomp.topolabel;
save(data_out_filename_ICA, 'icaunmixing', 'icatopolabel', 'eeg', 'eog_bp', 'browsed_data_raw', 'artifacts_data_raw');


