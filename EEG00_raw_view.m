function EEG00_raw_view(SBJ,view_previous, proc_id)
%% View raw data and mark epochs to toss
% INPUTS:
%   SBJ [str] - name of the subject to load
%   view_previous [0/1] - binary no/yes to load previous bad_epochs to view

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);

cfg=[];
cfg.dataset  = SBJ_vars.dirs.raw_filename;
cfg.demean   = 'yes';
% cfg.lpfilter = 'no';
cfg.hpfilter = 'yes';
cfg.hpfreq   = 0.5;
cfg.hpfiltord = 2;
%cfg.bpfilter = 'yes';
%cfg.bpfreq   = [0.5 40];%0.1 is too low for filter settings, 20 is too low to see muscle artifact, consider ditching filtering?
raw = ft_preprocessing(cfg);
%% Downsample
%Just looked over script a final time! This is one thing I wanted to check
%that its in the right spot that I forgot to mention -- I wanted it to go
%after filtering and trial cutting, but wasn't sure if this was too late!
%(April 13,2019)
if strcmp(proc_vars.resample_yn,'yes')
      cfg=[];
    cfg.resamplefs = proc_vars.resample_freq;
    raw = ft_resampledata(cfg, raw);
end
%% Plot PSDs for noise profile
fprintf('============== Plotting PSDs %s, %s ==============\n',SBJ,SBJ_vars.raw_file);
psd_dir = strcat(SBJ_vars.dirs.import,'raw_psds/');
if ~exist(psd_dir,'dir')
    mkdir(psd_dir);
end

fn_plot_PSD_1by1_save(raw.trial{1},raw.label,raw.fsample,...
    strcat(psd_dir,SBJ,'_raw_psd'),'png');

%% Plot and mark bad epochs
% Load cfg with plotting parameters
load([root_dir 'PRJ_Error_eeg/scripts/utils/cfg_plot_eeg.mat']);

% Load previous if available
out_fname = [SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat'];
if exist(out_fname)
    prev_exists = 1;
else
    prev_exists = 0;
end
if view_previous
    if prev_exists
        load(out_fname);
        cfg_plot.artfctdef.visual.artifact = bad_epochs;
    else
        warning('Previous bad_epochs from raw data cleaning doesnt exist! Please mark bad epochs.');
    end
end

browsed_raw = ft_databrowser(cfg_plot, raw);
bad_epochs  = browsed_raw.artfctdef.visual.artifact;

% prevent ft_databrowser bugs
fake_ep     = find(diff(bad_epochs,1,2)<10);
bad_epochs(fake_ep,:) = [];

%% Save out bad epochs
if prev_exists
    % Prevent overwriting most recent bad_epochs file
    mv_cmd = ['mv ' out_fname ' ' out_fname(1:end-4) '_bck.mat'];
    system(mv_cmd);
end
save(out_fname, 'bad_epochs');

end
