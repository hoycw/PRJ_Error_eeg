function SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_id,stat_id)
%% Run Mixed-Effects Linear model on ERPs from all SBJ and trials
%   NO RL MODEL! Only running with SBJ factor (random intercepts) to get
%   baseline model performance (R2 and AIC)
%   Select trials for conditions of interest
%       Optional: Find ERP, beta, or manual peak time to center stat window
%   Load single-trial ERP data and design matrix (only SBJ factor)
%   Run linear mixed effects model per time point or per electrode
%   Correct for multiple comparisons
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class, one cell per time point or channel
%   qvals [float array] - [n_regressors, n_chan/n_time] p values adjusted for multiple comparisons 
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)
%   time_vec [float array] - time points for each model run (same length as lme)
%   ch_list [cell array] - list of channels in analysis (for double checks)
%   reg_pk_time [float] - time used to center analysis window
%       NaN if run across time (st.measure = 'ts')

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data 
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.model_lab,'SBJonly'); error('Only model_lab = SBJonly for this script!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);

%% Load Behavior
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
for s = 1:numel(SBJs)
    % Load data
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    bhvs{s} = tmp.bhv;
    
    % Select Conditions of Interest
    full_cond_idx{s} = fn_condition_index(cond_lab, bhvs{s});
    bhv_fields = fieldnames(bhvs{s});
    orig_n_trials = numel(bhvs{s}.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhvs{s}.(bhv_fields{f_ix}))==orig_n_trials
            bhvs{s}.(bhv_fields{f_ix}) = bhvs{s}.(bhv_fields{f_ix})(full_cond_idx{s}~=0);
        end
    end
    n_trials(s) = numel(bhvs{s}.trl_n);
    
    clear tmp
end

%% Load Peak Timing Information
% This part of the analysis was not set up or run for SBJ_only!
% Can select peak time based on:
%   (1) positive or negative peak in window of ERP time series
%   (2) maximum model coefficient (absolute value) from previous model
%   (3) manually based on stat_vars
if strcmp(st.measure,'mean') && all(isfield(st,{'pk_reg_id','pk_stat_id','pk_an_id'}))
    error('not set up yet for null model');
    % Load previous stats
    tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.pk_stat_id '_' st.pk_an_id '.mat']);
    if numel(tmp.SBJs)~=numel(SBJs) || ~all(strcmp(tmp.SBJs,SBJs))
        error(['Not same SBJs in ' stat_id ' and ' st.pk_stat_id]);
    end
    
    % Obtain peak times for target regressor
    reg_ix = find(strcmp(reg_lab,st.pk_reg_id));
    pk_ts = nan(size(tmp.time_vec));
    for t_ix = 1:numel(tmp.time_vec)
        pk_ts(t_ix) = tmp.lme{t_ix}.Coefficients.Estimate(reg_ix+1);
    end
    [~,pk_ix] = max(abs(pk_ts));
    reg_pk_time = tmp.time_vec(pk_ix);
    
    % Adjust stat window based on peak time
    st.stat_lim = st.stat_lim+reg_pk_time;
else
    reg_pk_time = nan;
end

%% Load Data and Build Model
cfgs  = []; cfgs.latency = st.stat_lim;
sbj_factor  = zeros([sum(n_trials) 1]);
full_trl_ix = 0;
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',an_id,'.mat'],'roi');
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_roi = ft_selectdata(cfgs, roi);
    
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = st_roi.time{1};
        ch_list  = st_roi.label;
        if strcmp(st.measure,'ts')
            if numel(ch_list)>1; error('Why run st.measure = ts with more than one channel?'); end
            data  = nan([sum(n_trials) numel(time_vec)]);
        elseif strcmp(st.measure,'mean')
            data  = nan([sum(n_trials) numel(ch_list)]);
        else; error(['unknown st.measure: ' st.measure]);
        end
        
        % Create index for this SBJ in group design/data matrix
        sbj_idx = 1:n_trials(s);
    else
        sbj_idx = sum(n_trials(1:s-1))+1:sum(n_trials(1:s));
    end
    
    % Load and add data
    for trl_ix = 1:numel(st_roi.trial)
        full_trl_ix = full_trl_ix + 1;
        if strcmp(st.measure,'ts')
            data(full_trl_ix,:) = st_roi.trial{trl_ix};
        elseif strcmp(st.measure,'mean')
            data(full_trl_ix,:) = mean(st_roi.trial{trl_ix},2);
        else; error(['unknown st.measure: ' st.measure]);
        end
    end
    
    % Track SBJ in design matrix
    sbj_factor(sbj_idx) = s*ones([n_trials(s) 1]);
    
    clear roi st_roi
end

%% Run Linear Mixed Effects Model Over Time or Channels
fprintf('========================== Running Stats ==========================\n');
tic

% Build Model Table
tbl = table;
tbl.SBJ = categorical(sbj_factor);

% Create Model Formula
formula = 'ERP ~ 1 + (1|SBJ)';  % random intercepts for SBJ

% Run Model
if strcmp(st.measure,'ts')
    % Single channel over time
    lme = cell(size(time_vec));
    pvals = nan([numel(reg_lab) numel(time_vec)]);
    for t_ix = 1:numel(time_vec)
        tbl.ERP = data(:,t_ix);
        lme{t_ix} = fitlme(tbl,formula);
        pvals(:,t_ix) = lme{t_ix}.Coefficients.pValue;
    end
elseif strcmp(st.measure,'mean')
    % Single time point (mean in window) over channels
    lme = cell(size(ch_list));
    pvals = nan([numel(reg_lab) numel(ch_list)]);
    for ch_ix = 1:numel(ch_list)
        tbl.ERP = data(:,ch_ix);
        lme{ch_ix} = fitlme(tbl,formula);
        pvals(:,ch_ix) = lme{ch_ix}.Coefficients.pValue; % Take SBJ intercept
    end
else
    error(['Unknown st.measure: ' st.measure]);
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

fprintf('\t\t Stats Complete:');
toc

%% Save Results
stat_out_fname = [root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals','SBJs','time_vec','ch_list','reg_pk_time');

end
