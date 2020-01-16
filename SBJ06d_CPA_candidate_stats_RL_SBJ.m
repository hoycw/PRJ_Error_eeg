function SBJ06d_CPA_candidate_stats_RL_SBJ(SBJ,eeg_proc_id,cpa_id,an_id,stat_id)
%% Plot ERPs for single SBJ
% INPUTS:
%   conditions [str] - group of condition labels to segregate trials

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
cpa_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' cpa_id '_vars.m'];
eval(cpa_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Load data
load([SBJ_vars.dirs.proc SBJ '_' cpa_id '_' an_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' eeg_proc_id '_final.mat']);

model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Select Data to Model
full_cond_idx = fn_condition_index(cond_lab, bhv);
bhv_fields = fieldnames(bhv);
orig_n_trials = numel(bhv.trl_n);
for f_ix = 1:numel(bhv_fields)
    if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
        bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
    end
end

% Select time and trials in data
cfgs  = [];
cfgs.latency = st.stat_lim;
cfgs.trials  = find(full_cond_idx);
st_roi = ft_selectdata(cfgs, roi);
time_vec = st_roi.time{1};

% Create data matrix
if strcmp(st.measure,'ts')
    data  = nan([numel(bhv.trl_n) numel(time_vec)]);
elseif strcmp(st.measure,'mean')
    data  = nan(size(bhv.trl_n));
else; error(['unknown st.measure: ' st.measure]);
end
for trl_ix = 1:numel(st_roi.trial)
    if strcmp(st.measure,'ts')
        data(trl_ix,:) = st_roi.trial{trl_ix};
    elseif strcmp(st.measure,'mean')
        data(trl_ix) = mean(st_roi.trial{trl_ix});
    end
end

%% Build Model
% Load RL Model
load([SBJ_vars.dirs.proc SBJ '_model_' model_id '.mat']);

% Z-score SBJ model regressors
if st.z_reg
    for reg_ix = 1:numel(reg_lab)
        model(:,reg_ix) = ...
            (model(:,reg_ix)-nanmean(model(:,reg_ix)))./nanstd(model(:,reg_ix));
    end
end

%% Build Tabel and Run Model
fprintf('========================== Running Stats ==========================\n');
tic
% Build Model Table
tbl = table;
for reg_ix = 1:numel(reg_lab)
    %     if st.z_reg
    %         tbl.(reg_lab{reg_ix}) = ...
    %             (model(:,reg_ix)-nanmean(model(:,reg_ix)))./nanstd(model(:,reg_ix));
    %     else
    tbl.(reg_lab{reg_ix}) = model(:,reg_ix);
    %     end
end

% Create Model Formula
reg_formula = strjoin(reg_lab,' + ');
formula = ['ERP ~ ' reg_formula];

% Run Model
if strcmp(st.measure,'ts')
    lme = cell(size(time_vec));
    pvals = nan([numel(reg_lab) numel(time_vec)]);
    for t_ix = 1:numel(time_vec)
        tbl.ERP = data(:,t_ix);
        lme{t_ix} = fitlme(tbl,formula);
        pvals(:,t_ix) = lme{t_ix}.Coefficients.pValue(2:end);
    end
    
    % Correct for Multiple Comparisons
    if strcmp(st.mcp_method,'FDR')
        [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
        qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);
    else
        error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
    end
elseif strcmp(st.measure,'mean')
    lme = {};
    tbl.ERP = data;
    lme{1} = fitlme(tbl,formula);
    % No correction for multiple comparisons
    qvals = lme{1}.Coefficients.pValue(2:end);
end

fprintf('\t\t Stats Complete:');
toc

%% Save Results
stat_out_fname = [SBJ_vars.dirs.proc SBJ '_' stat_id '_' cpa_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals');

end
