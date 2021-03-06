function SBJ07c_CPA_candidate_ERP_GRP_stats_LME_RL(SBJ_id,eeg_proc_id,cpa_id,an_id,stat_id)
% Run Linear Mixed-Effects + Reinforcement Learning model on all SBJ candidate ICs
%   Only one channel for time series, multiple for topo
%   Can pull other CPA stats for peak timing

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Analysis and Plotting Parameters
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select Conditions of Interest
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]  = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Behavior
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
for s = 1:numel(SBJs)
    % Load behavior
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' eeg_proc_id '_final.mat'],'bhv');
    
    % Select Conditions of Interest
    full_cond_idx{s} = fn_condition_index(cond_lab, bhv);
    bhv_fields = fieldnames(bhv);
    orig_n_trials = numel(bhv.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
            bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx{s}~=0);
        end
    end
    n_trials(s) = numel(bhv.trl_n);
        
    clear bhv
end

%% Load Peak Timing Information
if strcmp(st.measure,'mean') && all(isfield(st,{'pk_reg_id','pk_stat_id','pk_an_id'}))
    % Load previous stats
    tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.pk_stat_id '_' cpa_id '_' st.pk_an_id '.mat']);
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
    st.stat_lim = st.stat_lim+reg_pk_time;
else
    reg_pk_time = nan;
end

%% Load Candidate IC, Build Model
cfgs  = []; cfgs.latency = st.stat_lim;
model = zeros([sum(n_trials) numel(reg_lab)]);
sbj_factor  = zeros([sum(n_trials) 1]);
full_trl_ix = 0;
for s = 1:numel(SBJs)
    % Load Candidate IC
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' ...
        SBJs{s} '_' cpa_id '_' an_id '.mat'],'roi');
    
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

        sbj_idx = 1:n_trials(s);
    else
        sbj_idx = sum(n_trials(1:s-1))+1:sum(n_trials(1:s));
    end
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' ...
        SBJs{s} '_model_' model_id '.mat']);
    
    % Z-score SBJ model regressors
    sbj_model = NaN(size(tmp.model));
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            sbj_model(:,reg_ix) = ...
                (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
        end
    else
        sbj_model = model;
    end
    model(sbj_idx,:) = sbj_model;
    
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
    
    % Track SBJ
    sbj_factor(sbj_idx) = s*ones([n_trials(s) 1]);
    
    clear roi st_roi
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
formula = ['ERP ~ ' reg_formula ' + (1|SBJ)'];

% Run Model
if strcmp(st.measure,'ts')
    lme = cell(size(time_vec));
    pvals = nan([numel(reg_lab) numel(time_vec)]);
    for t_ix = 1:numel(time_vec)
        tbl.ERP = data(:,t_ix);        
        lme{t_ix} = fitlme(tbl,formula);
        pvals(:,t_ix) = lme{t_ix}.Coefficients.pValue(2:end);
    end
elseif strcmp(st.measure,'mean')
    lme = cell(size(ch_list));
    pvals = nan([numel(reg_lab) numel(ch_list)]);
    for ch_ix = 1:numel(ch_list)
        tbl.ERP = data(:,ch_ix);
        lme{ch_ix} = fitlme(tbl,formula);
        % No correction for multiple comparisons
        pvals(:,ch_ix) = lme{ch_ix}.Coefficients.pValue(2:end);
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
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(stat_out_dir,'dir')
    [~] = mkdir(stat_out_dir);
end
stat_out_fname = [stat_out_dir SBJ_id '_' stat_id '_' cpa_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals','SBJs','time_vec','ch_list','reg_pk_time');

end
