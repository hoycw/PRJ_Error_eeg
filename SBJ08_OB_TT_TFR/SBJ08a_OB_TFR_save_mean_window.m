function SBJ08a_OB_TFR_save_mean_window(SBJ_id,proc_id,feat_id)
%% Save amplitude of condition-averaged TFR power in oddball (OB) task
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load and average single-trial TFR within condition
%   Average power within time-frequency window and across trials
%   Save amplitudes
% INPUTS:
%   SBJ_id [str]  - ID of subject list for group
%   proc_id [str] - ID of oddball preprocessing pipeline
%   feat_id [str] - ID of the feature extraction parameters to use
%       ft.an_id  = 'TFR_Fz_S2t1_dm2t0_fl05t20' || 'TFR_Pz_S2t1_dm2t0_fl05t20'
%       ft.grp_id = conditions to extract features (likely 'Tar','Odd','rare')
%       ft.measure:
%           'tfWin'- Mean window amplitude based on manual time-frequency range
% OUTPUTS:
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)
%   tfr_amp [float] - power amplitude in feature window

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data 
if ~contains(proc_id,'odd'); error('proc_id must be for oddball task!'); end
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' feat_id '_vars.m'];
eval(feat_vars_cmd);
if ~strcmp(ft.measure,'tfWin'); error('Must use TFR window for this script!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(ft.grp_id);

%% Load Data and Compute ERPs
tfr_amp = nan([numel(SBJs) numel(ft.name)]);
cfgs = [];
cfgs.avgovertime = 'yes';
cfgs.avgoverfreq = 'yes';
cfgs.avgoverrpt  = 'yes';
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_' proc_id '_' ft.an_id '.mat'],'tfr');
    
    for ft_ix = 1:numel(ft.name)
        % Select time and trials of interest
        cond_ix = find(strcmp(cond_lab,ft.cond{ft_ix}));
        cfgs.channel     = ft.chan{ft_ix};
        cfgs.latency     = ft.win_lim(ft_ix,:);
        cfgs.frequency   = ft.freq_lim(ft_ix,:);
        cfgs.trials      = find(fn_condition_index(cond_lab, bhv)==cond_ix);
        
        ft_tfr = ft_selectdata(cfgs, tfr);
        tfr_amp(s,ft_ix) = ft_tfr.powspctrm;
    end
    
    clear bhv tfr ft_tfr
end

%% Save Results
feat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(feat_out_dir,'dir')
    mkdir(feat_out_dir);
end
stat_out_fname = [feat_out_dir SBJ_id '_' feat_id '_' proc_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','SBJs','tfr_amp');

end
