function EEG00_raw_view(SBJ,view_previous)
%% View raw data and mark epochs to toss
% INPUTS:
%   SBJ [str] - name of the subject to load
%   view_previous [0/1] - binary no/yes to load previous bad_epochs to view

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
cfg.bpfreq   = [0.5 40];%0.1 is too low for filter settings, 20 is too low to see muscle artifact, consider ditching filtering?
raw = ft_preprocessing(cfg);

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
bad_epochs = browsed_raw.artfctdef.visual.artifact;

%% Save out bad epochs
if prev_exists
    % Prevent overwriting most recent bad_epochs file
    mv_cmd = ['mv ' out_fname ' ' out_fname(1:end-4) '_bck.mat'];
    system(mv_cmd);
end
save(out_fname, 'bad_epochs');

end
