function SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_id,stat_id)
%% Run Mixed-Effects Linear model on singel-trial EEG from all SBJ and trials
%   Either a single channel across time or single time point across
%   channels (e.g., average in time window)
% COMPUTATIONS:
%   Select trials for conditions of interest
%       Optional: Find ERP, beta, or manual peak time to center stat window
%   Load single-trial ERP data and design matrix (model regressors, SBJ factor)
%       Optional: z-score model regressors within SBJ
%   Compute and plot group concatenated model (design matrix) and correlations
%   Run linear mixed effects model per time point or per electrode
%   Correct for multiple comparisons (FDR for regressors, time points or channels)
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
if strcmp(st.measure,'mean')
    if all(isfield(st,{'pk_stat_cond','pk_erp_cond','pk_lim','pk_sign'}))
        % (1) Find ERP peak
        % Load ERP
        tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.pk_stat_cond '_' st.pk_an_id '.mat']);
        
        % Select Time Windows within Conditions
        [pk_cond_lab] = fn_condition_label_styles(st.pk_stat_cond);
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
        % Check these stats were run on the same subjects
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
model = zeros([sum(n_trials) numel(reg_lab)]);
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
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' st.model_id '.mat']);
    
    % Z-score SBJ model regressors
    sbj_model = NaN([sum(full_cond_idx{s}~=0) size(tmp.model,2)]);
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            if strcmp(st.model_cond,'DifFB') && ~strcmp(st.stat_cond,'DifFB')
                % Model was run on all trials, but using subset here
                sbj_model(:,reg_ix) = ...
                    (tmp.model(full_cond_idx{s}~=0,reg_ix)-nanmean(tmp.model(full_cond_idx{s}~=0,reg_ix)))./...
                    nanstd(tmp.model(full_cond_idx{s}~=0,reg_ix));
            elseif strcmp(st.model_cond,st.stat_cond)
                % Model and Stats on same trials, order matches
                sbj_model(:,reg_ix) = ...
                    (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./...
                    nanstd(tmp.model(:,reg_ix));
            else
                error(['Mismatch between st.model_cond ' st.model_cond ' and st.stat_cond ' st.stat_cond]);
            end
        end
    else
        sbj_model = tmp.model(full_cond_idx{s}~=0,:);
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
    
    % Track SBJ in design matrix
    sbj_factor(sbj_idx) = s*ones([n_trials(s) 1]);
    
    clear roi st_roi
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
fig_dir = [stat_out_dir st.model_id '_' st.stat_cond '_plots/'];
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

%% Run Linear Mixed Effects Model Over Time or Channels
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
    % Single channel over time
    lme = cell(size(time_vec));
    pvals = nan([numel(reg_lab) numel(time_vec)]);
    for t_ix = 1:numel(time_vec)
        tbl.ERP = data(:,t_ix);
        lme{t_ix} = fitlme(tbl,formula);
        pvals(:,t_ix) = lme{t_ix}.Coefficients.pValue(2:end);
    end
elseif strcmp(st.measure,'mean')
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
