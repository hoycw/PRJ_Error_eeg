function SBJ04c_ERP_grp_stats_LME_mean_window(SBJ_id,proc_id,an_id,stat_id)
%% Run Mixed-Effects Linear model on condition-averaged ERPs from all SBJs
%   Averages ERPs (condition averaged data per SBJ) in time window and across channels
% COMPUTATIONS:
%   Select trials for conditions of interest
%       Optional: Find ERP, beta, or manual peak time to center stat window
%   Load ERP data and design matrix (model regressors, SBJ factor)
%       Optional: z-score model regressors within SBJ
%   Compute and plot group concatenated model (design matrix) and correlations
%   Run linear mixed effects model per time point or per electrode
%   Correct for multiple comparisons (FDR for regressors)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
%       st.measure = 'erp_mean' averages across single SBJ ERPs
%           this is default because it best approximates previous literature
%                  = 'mean' (depricated) averages across single trial amplitudes 
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class, one cell per channel
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
if ~strcmp(st.measure,'erp_mean') || ~strcmp(st.an_style,'lme'); error('This LME should only be run on mean window!'); end

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
% Can select peak time based on:
%   (1) positive or negative peak in window of ERP time series
%   (2) maximum model coefficient (absolute value) from previous model
%   (3) manually based on stat_vars
if any(strcmp(st.measure,{'mean','erp_mean'}))
    if all(isfield(st,{'pk_cond_grp','pk_erp_cond','pk_lim','pk_sign'}))
        % (1) Find ERP peak
        % Load ERP
        tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.pk_cond_grp '_' st.pk_an_id '.mat']);
        
        % Select Time Windows within Conditions
        [pk_cond_lab] = fn_condition_label_styles(st.pk_cond_grp);
        cfgs = []; cfgs.latency = st.pk_lim;
        st_erp = ft_selectdata(cfgs,tmp.er_grp{strcmp(pk_cond_lab,st.pk_erp_cond)});
        
        % Obtain peak times in window
        pk_ts = st_erp.avg;
        [~,pk_ix] = max(pk_ts*st.pk_sign);
        reg_pk_time = st_erp.time(pk_ix);
    elseif all(isfield(st,{'pk_reg_id','pk_stat_id'}))
        % (2) Find regression beta peak
        % Load previous stats
        tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.pk_stat_id '_' st.pk_an_id '.mat']);
        if numel(tmp.SBJs)~=numel(SBJs) || ~all(strcmp(tmp.SBJs,SBJs))
            error(['Not same SBJs in ' stat_id ' and ' st.pk_stat_id]);
        end
        
        % Obtain peak times for target regressor
        real_st = st;
        eval(['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' real_st.pk_stat_id '_vars.m']);
        pk_st = st; st = real_st;
        [pk_reg_lab, ~, ~, ~] = fn_regressor_label_styles(pk_st.model_lab);
        reg_ix = find(strcmp(pk_reg_lab,st.pk_reg_id));
        pk_ts = nan(size(tmp.time_vec));
        for t_ix = 1:numel(tmp.time_vec)
            pk_ts(t_ix) = tmp.lme{t_ix}.Coefficients.Estimate(reg_ix+1);
        end
        [~,pk_ix] = max(abs(pk_ts));
        reg_pk_time = tmp.time_vec(pk_ix);
    elseif isfield(st,'pk_center')
        % (3) Use manual peak setting
        reg_pk_time = st.pk_center;
    end
    
    % Adjust stat window based on peak time
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
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' st.model_id '.mat']);
    
    % Z-score SBJ model regressors
    sbj_model = NaN([sum(full_cond_idx{s}~=0) size(tmp.model,2)]);
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            sbj_model(:,reg_ix) = ...
                (tmp.model(full_cond_idx{s}~=0,reg_ix)-nanmean(tmp.model(full_cond_idx{s}~=0,reg_ix)))./...
                nanstd(tmp.model(full_cond_idx{s}~=0,reg_ix));
        end
    else
        sbj_model = tmp.model(full_cond_idx{s}~=0,:);
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
    
    % Track SBJ in design matrix
    sbj_factor = [sbj_factor; s*ones([numel(cond_lab) 1])];
    
    clear tmp roi st_roi sbj_model trials cond_idx cond_trial_ix
end

%% Compute and plot correlations between regressors
reg_corr = corr(model,'rows','complete');
rvals = reshape(triu(reg_corr,1),[numel(reg_corr) 1]);
[~,max_r_ix] = max(abs(rvals));
if st.z_reg
    reg_vifs = fn_variance_inflation_factor(model);
else
    reg_vifs = fn_variance_inflation_factor(zscore(model));
end
[max_vif,max_vif_ix] = max(reg_vifs);

% Create figure directory
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
fig_dir = [stat_out_dir st.model_id '_' st.stat_cond '_erpMW_plots/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Plot design matrix
fig_name = [SBJ_id '_' st.model_id '_' st.stat_cond '_design'];
figure('Name',fig_name);
imagesc(model);
xticklabels(reg_lab);
colorbar;
saveas(gcf,[fig_dir fig_name '.png']);

% Plot regressor correlation matrix
fig_name = [SBJ_id '_' st.model_id '_' st.stat_cond '_design_corr_VIFs'];
figure('Name',fig_name);
subplot(1,2,1);
imagesc(reg_corr);
set(gca,'XLim',[0.5 numel(reg_lab)+0.5]);
set(gca,'XTick',1:numel(reg_lab));
set(gca,'XTickLabels',reg_lab);
set(gca,'YLim',[0.5 numel(reg_lab)+0.5]);
set(gca,'YTick',1:numel(reg_lab));
set(gca,'YTickLabels',reg_lab);
colorbar;
title(['Max r = ' num2str(rvals(max_r_ix),'%.3f')]);
set(gca,'FontSize',16);

% Plot VIFs
subplot(1,2,2);
bars = bar(1:numel(reg_lab),reg_vifs);
bars.FaceColor = 'flat';
bars.CData(2,:) = [.5 0 .5];
for reg_ix = 1:numel(reg_lab)
    if reg_vifs(reg_ix) >= 5 && reg_vifs(reg_ix) < 10
        bars.CData(reg_ix,:) = [1 .5 0];
    elseif reg_vifs(reg_ix) >= 10
        bars.CData(reg_ix,:) = [1 0 0];
    else
        bars.CData(reg_ix,:) = [0 0 0];
    end
end
set(gca,'XLim',[0.5 numel(reg_lab)+0.5]);
set(gca,'XTick',1:numel(reg_lab));
set(gca,'XTickLabels',reg_lab);
set(gca,'FontSize',16);
title(['Max VIF: ' reg_lab{max_vif_ix} '=' num2str(max_vif,'%.2f')]);
saveas(gcf,[fig_dir fig_name '.png']);

%% Run Linear Mixed Effects Model Over Channels
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
formula = ['ERP ~ ' reg_formula ' + (1|SBJ)'];  % random intercepts for SBJ

% Run Model
if strcmp(st.measure,'ts')
    error('must use st.measure = erp_mean for this script!');
elseif any(strcmp(st.measure,{'mean','erp_mean'}))
    % Single time point (mean in window) over channels
    lme = cell(size(ch_list));
    pvals = nan([numel(reg_lab) numel(ch_list)]);
    for ch_ix = 1:numel(ch_list)
        tbl.ERP = data(:,ch_ix);
        lme{ch_ix} = fitlme(tbl,formula);
        pvals(:,ch_ix) = lme{ch_ix}.Coefficients.pValue(2:end); % Skip intercept
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
