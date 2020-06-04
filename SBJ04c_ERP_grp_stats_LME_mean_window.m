function SBJ04c_ERP_grp_stats_LME_mean_window(SBJ_id,proc_id,an_id,stat_id)
% Run Mixed-Effects Linear model on condition-averaged SBJ ERPs
%   Only for one channel now...
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class for each time point

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
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
if ~strcmp(st.measure,'erp_mean') || ~strcmp(st.an_style,'lme'); error('This LME should only be run on mean window!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

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

%% Load Peak Timing Information
if any(strcmp(st.measure,{'mean','erp_mean'}))
    if all(isfield(st,{'pk_trial_cond','pk_erp_cond','pk_lim','pk_sign'}))    % ERP peak
        % Load ERP
        tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.pk_trial_cond '_' st.pk_an_id '.mat']);
        
        % Select Time Windows
        [pk_cond_lab] = fn_condition_label_styles(st.pk_trial_cond);
        cfgs = []; cfgs.latency = st.pk_lim;
        st_erp = ft_selectdata(cfgs,tmp.er_grp{strcmp(pk_cond_lab,st.pk_erp_cond)});
        
        % Obtain peak times in window
        pk_ts = st_erp.avg;
        [~,pk_ix] = max(pk_ts*st.pk_sign);
        reg_pk_time = st_erp.time(pk_ix);
    elseif all(isfield(st,{'pk_reg_id','pk_stat_id'}))   % regression beta peak
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
    elseif isfield(st,'pk_center')      % manual peak setting
        reg_pk_time = st.pk_center;
    end
    st.stat_lim = st.stat_lim+reg_pk_time;
else
    reg_pk_time = nan;
end

%% Load Data and Build Model
cfgs  = []; cfgs.latency = st.stat_lim;
model = zeros([numel(cond_lab)*numel(SBJs) numel(reg_lab)]);
data  = zeros([numel(cond_lab)*numel(SBJs) 1]);
sbj_factor  = [];
full_model_ix = 0;
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
        erps  = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
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
    
    % Compute ERPs and average model within condition
    cond_idx = fn_condition_index(cond_lab, bhvs{s});
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        % Compute ERP
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for trl_ix = 1:numel(cond_trial_ix)
            trials(:,trl_ix,:) = st_roi.trial{cond_trial_ix(trl_ix)};
        end
        
        % Average data and model
        full_model_ix = full_model_ix + 1;
        data(full_model_ix)     = mean(mean(trials,2));
        model(full_model_ix, :) = mean(sbj_model(cond_trial_ix,:),1);
    end
    
    % Track SBJ
    sbj_factor = [sbj_factor; s*ones([numel(cond_lab) 1])];
    
    clear tmp roi st_roi sbj_model trials cond_idx cond_trial_ix
end

%% Compute and plot correlations between regressors
reg_corr = corr(model,'rows','complete');

% Create figure directory
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
fig_dir = [stat_out_dir model_id '_erp_plots/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Plot design matrix
fig_name = [SBJ_id '_' model_id '_design'];
figure('Name',fig_name);
imagesc(model);
xticklabels(reg_lab);
colorbar;
saveas(gcf,[fig_dir fig_name '.png']);

% Plot regressor correlation matrix
fig_name = [SBJ_id '_' model_id '_design_corr'];
figure('Name',fig_name);
imagesc(reg_corr);
xticklabels(reg_lab);
yticklabels(reg_lab);
colorbar;
saveas(gcf,[fig_dir fig_name '.png']);

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
%     lme = cell(size(time_vec));
%     pvals = nan([numel(reg_lab) numel(time_vec)]);
%     for t_ix = 1:numel(time_vec)
%         tbl.ERP = data(:,t_ix);
%         lme{t_ix} = fitlme(tbl,formula);
%         pvals(:,t_ix) = lme{t_ix}.Coefficients.pValue(2:end);
%     end
elseif any(strcmp(st.measure,{'mean','erp_mean'}))
    lme = cell(size(ch_list));
    pvals = nan([numel(reg_lab) numel(ch_list)]);
    for ch_ix = 1:numel(ch_list)
        tbl.ERP = data(:,ch_ix);
        lme{ch_ix} = fitlme(tbl,formula);
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
stat_out_fname = [stat_out_dir SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals','SBJs','time_vec','ch_list','reg_pk_time');

end
