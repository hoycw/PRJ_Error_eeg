function EEG00_raw_view(SBJ)
%% Only for viewing raw data and marking epochs to toss

%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath(genpath([root_dir 'PRJ_Error_eeg/scripts/']));
addpath(ft_dir);
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

cfg=[];
cfg.dataset  = SBJ_vars.dirs.raw_filename;
cfg.lpfilter = 'no';
cfg.hpfilter = 'no';
cfg.bpfilter = 'yes';
%low pass filter at 8 hz
cfg.bpfreq   = [0.1 200];
%cfg.trialdef.eventtype  = 'STATUS';
%cfg.trialdef.eventvalue = 2;
%cfg.trialdef.prestim    = -.3;
%cfg.trialdef.poststim   = 1;
%cfg.trialfun            = 'tt_trialfun';
%cfg = ft_definetrial(cfg);
dataraw = ft_preprocessing(cfg);
cfg.viewmode = 'vertical';
browsed_data_raw = ft_databrowser(cfg, dataraw);
%data = ft_rejectartifact(cfg, dataraw);
save(data_out_filename_rawdatainspect, 'browsed_data_raw', 'dataraw');

end