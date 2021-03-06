function SBJ00_raw_view(SBJ,view_previous, proc_id, block)
%% View raw data and mark epochs to exclude from ICA
% INPUTS:
%   SBJ [str] - name of the subject to load
%   view_previous [0/1] - binary no/yes to load previous bad_epochs to view
%   proc_id [str] - name of processing pipeline
%   block [int] - which block number to run
% OUTPUTS:
%   bad_epochs [Nx2 int] - [start stop] samples for N bad epochs
%       if previous bad_epochs exists, moves to bad_epochs_bck.mat

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

if numel(SBJ_vars.block_name)>1
        block_suffix = ['_' SBJ_vars.block_name{block}];
else
        block_suffix = '';
end
cfg=[];
cfg.dataset  = SBJ_vars.dirs.raw_filename{block};
cfg.demean   = 'yes';
cfg.hpfilter = 'yes';
cfg.hpfreq   = 0.5;
cfg.hpfiltord = 2;
raw = ft_preprocessing(cfg);

%% Downsample
if strcmp(proc.resample_yn,'yes')
    cfg = [];
    cfg.resamplefs = proc.resample_freq;
    raw = ft_resampledata(cfg, raw);
end

%% Plot PSDs for noise profile
fprintf('============== Plotting PSDs %s, %s ==============\n',SBJ,SBJ_vars.raw_file{block});
psd_dir = strcat(SBJ_vars.dirs.import,'raw_psds/');
if ~exist(psd_dir,'dir')
    mkdir(psd_dir);
end

fn_plot_PSD_1by1_save(raw.trial{1},raw.label,raw.fsample,...
    strcat(psd_dir,SBJ,'_raw_psd',block_suffix),'png');

%% Plot and mark bad epochs
% Load cfg with plotting parameters
load([root_dir 'PRJ_Error_eeg/scripts/utils/cfg_plot_eeg.mat']);

% Load previous if available
out_fname = [SBJ_vars.dirs.events SBJ '_raw_bad_epochs' block_suffix '.mat'];
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
% Prevent ft_databrowser tiny epoch bugs
bad_epochs(diff(bad_epochs,1,2)<10,:) = [];

%% Save out bad epochs
if prev_exists
    % Prevent overwriting most recent bad_epochs file
    mv_cmd = ['mv ' out_fname ' ' out_fname(1:end-4) '_bck.mat'];
    system(mv_cmd);
end
save(out_fname, 'bad_epochs');

end
