function SBJ04e_ERP_print_RL_model_comparison_win(SBJ_id,an_id,stat_ids,null_id,mean_win,varargin)
%% Print AIC model performance for different models averaged in windows
%   Also adds relatively likelihoods in legend
%   Only for single channel right now...
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   an_id [str] - ID of the analysis parameters to use
%   stat_ids [cell array] - string IDs of the stats parameters to compare
%       Try to add the main model last, since it will overlay others
%       3 model colors: magenta, lime green, black (otherwise R,G,B,etc.)
%   null_id [str] - ID of the SBJonly baseline model to compare
%   mean_win [2x1 array] - [start, stop] times in sec for averaging AIC and computing relative likelihoods
%   varargin:
%       rm_null [0/1] - binary flag to subtract out null model with only random intercepts
%           default: 0
%       mean_reg [cell array] - {reg_lab, stat_id} to center 50ms win for mean AIC and relative likelihoods
% OUTPUTS:
%   prints AIC values avearged in windows

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'mean_reg')
            mean_reg = varargin{v+1};
        elseif strcmp(varargin{v},'rm_null')
            rm_null = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('rm_null','var'); rm_null = 0; end

%% Analysis and Plotting Parameters
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Load stat parameters and check compatibility
sts = cell(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_ids{st_ix} '_vars.m'];
    eval(stat_vars_cmd);
    sts{st_ix} = st;
    if ~strcmp(st.measure,'ts'); error('this script is for time series!');end
    
    % Check alignment of time windows and measurements
    if st_ix>1
        if any(sts{1}.stat_lim ~= sts{st_ix}.stat_lim)
            error('st.stat_lim not aligned!');
        end
        if ~strcmp(sts{1}.measure, sts{st_ix}.measure)
            error('st.measure not the same!');
        end
    end
    clear st stat_vars_cmd
end

% Load SBJonly null model
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' null_id '_vars.m'];
eval(stat_vars_cmd);
null_st = st;

% Check compatibility of null model
if any(sts{1}.stat_lim ~= null_st.stat_lim)
    error('null_st.stat_lim not aligned!');
end
if ~strcmp(sts{1}.measure, null_st.measure)
    error('null_st.measure not the same!');
end
clear st stat_vars_cmd

% Get Plotting Parameters
load([root_dir 'PRJ_Error_EEG/data/' SBJs{1} '/04_proc/' SBJs{1} '_' an_id '.mat'],'roi');
cfgs = []; cfgs.latency = sts{1}.stat_lim;
st_roi = ft_selectdata(cfgs, roi);
st_time_vec = st_roi.time{1};
ch_list = st_roi.label;
if exist('mean_win','var')
    [~, start_ix] = min(abs(st_time_vec-mean_win(1)));
    [~, stop_ix] = min(abs(st_time_vec-mean_win(2)));
    aic_win_idx = start_ix:stop_ix;
else
    aic_win_idx = 1:numel(st_time_vec);
end

%% Load Models
% Load real models
lmes = cell([numel(stat_ids) numel(st_time_vec)]);
for st_ix = 1:numel(stat_ids)
    tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_ids{st_ix} '_' an_id '.mat']);
    lmes(st_ix,:) = tmp.lme;
end

% Load null model
null_aic = zeros(size(st_time_vec));
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' null_id '_' an_id '.mat']);
for t_ix = 1:numel(st_time_vec)
    null_aic(t_ix) = tmp.lme{t_ix}.ModelCriterion.AIC;
end

% Grab peak time
if exist('mean_reg','var')
    % Select regressor and stat_id
    stat_id_ix = strcmp(stat_ids,mean_reg{2});
    [reg_lab, ~, ~, ~]  = fn_regressor_label_styles(sts{stat_id_ix}.model_lab);
    reg_ix = find(strcmp(reg_lab,mean_reg{1}));
    % Select betas
    beta_ts = nan(size(st_time_vec));
    for t_ix = 1:numel(st_time_vec)
        beta_ts(t_ix) = lmes{stat_id_ix,t_ix}.Coefficients.Estimate(reg_ix+1);
    end
    % Obtain peak times for target regressor
    [~,pk_ix] = max(abs(beta_ts));
    [~,start_ix] = min(abs(st_time_vec-st_time_vec(pk_ix)+0.025));
    [~,stop_ix]  = min(abs(st_time_vec-st_time_vec(pk_ix)-0.025));
    aic_win_idx = start_ix:stop_ix;
    mean_win = [st_time_vec(start_ix) st_time_vec(stop_ix)];
end

%% Compute model performance data
% Collect AIC
aics = NaN([numel(stat_ids) numel(st_time_vec)]);
mean_aic = NaN(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    for t_ix = 1:numel(st_time_vec)
        aics(st_ix,t_ix) = lmes{st_ix,t_ix}.ModelCriterion.AIC;
    end
    if rm_null
        aics(st_ix,:) = aics(st_ix,:) - null_aic;
    end
    mean_aic(st_ix) = nanmean(aics(st_ix,aic_win_idx));
end

mean_aic_null = nanmean(null_aic(aic_win_idx));

%% Compute relative likelihoods and pirnt AIC
fprintf('%s mean window: [%.3f to %.3f]\n',an_id,mean_win(1),mean_win(2));
aic_min = min(mean_aic);
rel_lik = nan(size(stat_ids));
st_leg  = cell(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    rel_lik(st_ix) = exp((aic_min-mean_aic(st_ix))/2);
    fprintf('\t%s\n',[stat_ids{st_ix} ' (mean=' num2str(round(mean_aic(st_ix)))...
        '; RL=' num2str(rel_lik(st_ix),'%.2f') ')']);
end

% Compute for null model
rel_lik_null  = exp((aic_min-mean_aic_null)/2);

fprintf(2,'\t%s\n',[null_id ' (mean=' num2str(round(mean_aic_null)) ...
    '; RL=' num2str(rel_lik_null,'%.2f') ')']);


end
