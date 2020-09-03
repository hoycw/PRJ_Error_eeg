function SBJ05d_TFR_grp_stats_LME_RL(SBJ_id,proc_id,an_id,stat_id)
%% Run Linear Mixed-Effects multiple regression model on single-trial
% time-frequency power for all SBJ and trials
%   Only for single channels
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load single-trial TFR power data and design matrix (model regressors, SBJ factor)
%       Optional: z-score model regressors within SBJ
%   Run LME RL model per time-frequency point
%   Correct for multiple comparisons (FDR for regressors, times, frequencies)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class for each [frequency, time] point
%   qvals [float array] - [n_regressors, n_freq, n_time] p values adjusted for multiple comparisons 
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)

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
if ~strcmp(st.an_style,'lme'); error('stat_id not using lme!'); end
if ~strcmp(st.measure,'ts'); error('why are you trying to average across TFR?'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

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

%% Load Data and Build Model
cfgs  = []; cfgs.latency = st.stat_lim;
model = zeros([sum(n_trials) numel(reg_lab)]);
sbj_factor  = zeros([sum(n_trials) 1]);
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_' proc_id '_' an_id '.mat'],'tfr');
    if numel(tfr.label)>1; error('assuming single channel for now!'); end
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_tfr = ft_selectdata(cfgs, tfr);
    
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = st_tfr.time;
        fois     = st_tfr.freq;
        if strcmp(st.measure,'ts')
            data  = nan([sum(n_trials) numel(fois) numel(time_vec)]);
        elseif strcmp(st.measure,'mean')
            error('why are you trying to average across TFR?');
%             data  = nan([sum(n_trials) numel(fois)]);
        else; error(['unknown st.measure: ' st.measure]);
        end
        
        sbj_idx = 1:n_trials(s);
    else
        sbj_idx = sum(n_trials(1:s-1))+1:sum(n_trials(1:s));
    end
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
    
    % Z-score SBJ model regressors
    sbj_model = NaN(size(tmp.model));
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            sbj_model(:,reg_ix) = ...
                (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
        end
    else
        sbj_model = tmp.model;
    end
    model(sbj_idx,:) = sbj_model;
    
    % Load and add power data
    if strcmp(st.measure,'ts')
        data(sbj_idx,:,:) = squeeze(st_tfr.powspctrm);
    elseif strcmp(st.measure,'mean')
%         data(sbj_idx,:,:) = squeeze(mean(st_tfr.powspctrm,4));
    else; error(['unknown st.measure: ' st.measure]);
    end
    
    % Track SBJ in design matrix
    sbj_factor(sbj_idx) = s*ones([n_trials(s) 1]);
    
    clear tfr st_tfr
end

%% Build table
fprintf('========================== Running Stats ==========================\n');
tic
% Build Model Table
tbl = table;
for reg_ix = 1:numel(reg_lab)
    tbl.(reg_lab{reg_ix}) = model(:,reg_ix);
end
tbl.SBJ = categorical(sbj_factor);

% Create Model Formula
reg_formula = strjoin(reg_lab,' + ');
formula = ['TFR ~ ' reg_formula ' + (1|SBJ)'];  % random intercepts for SBJ

% Run Model for each time and frequency
if strcmp(st.measure,'ts')
    lme = cell([numel(fois) numel(time_vec)]);
    pvals = nan([numel(reg_lab) numel(fois) numel(time_vec)]);
    for f_ix = 1:numel(fois)
        for t_ix = 1:numel(time_vec)
            tbl.TFR = data(:,f_ix,t_ix);
            lme{f_ix,t_ix} = fitlme(tbl,formula);
            pvals(:,f_ix,t_ix) = lme{f_ix,t_ix}.Coefficients.pValue(2:end);
        end
    end
    
    % Correct for Multiple Comparisons (regressors, times, frequencies)
    if strcmp(st.mcp_method,'FDR')
        [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2)*size(pvals,3) 1]));
        qvals = reshape(qvals,[size(pvals,1) size(pvals,2) size(pvals,3)]);
    else
        error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
    end
elseif strcmp(st.measure,'mean')
%     lme =cell(size(fois));
%     for f_ix = 1:numel(fois)
%         tbl.TFR = data(:,f_ix);
%         lme{f_ix} = fitlme(tbl,formula);
%         pvals(:,f_ix) = lme{f_ix}.Coefficients.pValue(2:end);
%     end
%     
%     % Correct for Multiple Comparisons
%     if strcmp(st.mcp_method,'FDR')
%         [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
%         qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);
%     else
%         error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
%     end
end

fprintf('\t\t Stats Complete:');
toc

%% Save Results
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(stat_out_dir,'dir')
    [~] = mkdir(stat_out_dir);
end
stat_out_fname = [stat_out_dir SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals','SBJs');

end
