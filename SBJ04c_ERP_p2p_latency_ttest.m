function SBJ04c_ERP_p2p_latency_ttest(SBJ_id,an_id,pk_stat_id,SBJ_norm)
%% Run paired samples t-test on FRN peak latency
%   Meant only for testing easy/hard surprise latencies
%   Must run SBJ04c_ERP_grp_stats_LME_P2P first to obtain peak data
%   Only for single channel
% COMPUTATIONS:
%   Load single-trial design matrix (model regressors) and average within condtion
%   Load peak times identified in peak-to-peak FRN LME analysis
%   LME multiple regression predicting peak latency using model regressors
%   Correct for multiple comparisons (FDR for regressors)
%   Scatter plot of latencies with simple linear fit for visualization
%       Plotting Option: normalize latencies within SBJ by subtracting mean latency across conditions
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   pk_stat_id [str] - ID of the peak-to-peak stats analysis to provide peak data
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure (and prints outcome statistics)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Analysis and Plotting Parameters
if ~contains(pk_stat_id,'_EHSu_lme_p2pFRN')
    error('This script is only meant for peak latency t-test between easy/hard surprise');
end
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' pk_stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.measure,'p2p') || ~strcmp(st.an_style,'lme'); error('run only for p2p LME analyses!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);

%% Load Peak Times
% Can't load P2P results from only easy/hard surprise, so select
%   conditions after loading full DifFB results
pk_times_id = strrep(pk_stat_id,'EHSu','DifFB');

% Load P2P results (peak times)
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' pk_times_id '_' an_id '.mat']);
if numel(ch_list)>1; error('only plotting for 1 channel in this script!'); end
if st.pk_sign(2)~=-1; error('second peak is not negative!'); end

% Take second (negative) peak; assume 1 channel
pk_data = squeeze(pk_times(:,:,1,2));

% Normalize peak latencies within SBJ
if SBJ_norm
    pk_data = pk_data-mean(pk_data,1);
    norm_str = '_SBJnorm';
else
    norm_str = '';
end

% Subselect conditions after normalizing across all
[gen_cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.model_cond);
cond_idx = contains(gen_cond_lab,cond_lab);
pk_data = pk_data(cond_idx,:);

%% Compute paired samples t-test for latency across conditions
[~,pval,~,stats] = ttest(pk_data(1,:),pk_data(2,:));

% Print results
fprintf('%s (n=%d) %s mean +/- STD:\n\t%s = %f +/- %f\n\t%s = %f +/- %f\n',...
    SBJ_id, numel(SBJs), norm_str, cond_lab{1}, nanmean(pk_data(1,:)), nanstd(pk_data(1,:)),...
    cond_lab{2}, nanmean(pk_data(2,:)), nanstd(pk_data(2,:)));
fprintf('p = %f, t=%f, df=%d\n',pval, stats.tstat, stats.df);

end
