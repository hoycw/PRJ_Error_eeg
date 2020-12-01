function SBJ06c_OB_ERP_TT_grp_stats_reg_mean_window(SBJ_id,tt_proc_id,ob_proc_id,an_id,stat_id,varargin)
%% Run multiple regression using OB ERP features to predict TT mean window amplitude
%   Oddball ERPs only yield one feature per SBJ, so can only predict single TT condition/average
%   Averages ERPs (averaged within condition, then across conditions per SBJ) in time window
% COMPUTATIONS:
%   Select trials for conditions of interest
%       Optional: Find ERP, beta, or manual peak time to center stat window
%   Load ERP data and predictor matrix (OB ERP features, SBJ factor)
%       Optional: z-score regressors within SBJ
%   Compute and plot group concatenated model (design matrix) and correlations
%   Run linear mixed effects model per time point or per electrode
%   Correct for multiple comparisons (FDR for regressors)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   tt_proc_id [str] - ID of target time preprocessing pipeline
%   ob_proc_id [str] - ID of oddball preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
%       st.model = oddball feat_id for OB ERP features
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

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');    fig_vis = 'on'; end
if ~exist('fig_ftype','var');  fig_ftype = 'png'; end
if ~exist('save_fig','var');   save_fig = 1; end

%% Load Data 
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.measure,'erp_mean') || ~strcmp(st.an_style,'reg'); error('This regression should only be run on mean window!'); end
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' st.model_lab '_vars.m'];
eval(feat_vars_cmd);
if ~any(strcmp(ft.grp_id,{'rare','Odd','Tar'})); error('Features should be oddball conditions!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);

%% Load Behavior
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
for s = 1:numel(SBJs)
    % Load data
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' tt_proc_id '_final.mat'],'bhv');
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

%% Load Oddball ERP Features
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'SBJs');
if numel(SBJs)~=numel(tmp.SBJs) || ~all(strcmp(SBJs,tmp.SBJs)); error('SBJ mismatch!'); end
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'ft_amp','ft_times');

% Z-score feature predictors
model = ft_amp;
if st.z_reg
    model = zscore(model);
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
        model_id_strs = strsplit(st.pk_stat_id,'_');
        [reg_lab, ~, ~, ~]  = fn_regressor_label_styles(model_id_strs{1});
        
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
data  = zeros([numel(SBJs) 1]);
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
    
    % Compute ERPs and average model within condition
    cond_idx = fn_condition_index(cond_lab, bhvs{s});
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        % Compute ERP
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for trl_ix = 1:numel(cond_trial_ix)
            trials(:,trl_ix,:) = st_roi.trial{cond_trial_ix(trl_ix)};
        end
        
        % Average data
        data(s) = mean(mean(trials,2));
    end
        
    clear tmp roi st_roi trials cond_idx cond_trial_ix
end

%% Run Linear Multiple Regression
fprintf('========================== Running Stats ==========================\n');
% Build Model Table
tbl = table;
for ft_ix = 1:numel(ft.name)
    tbl.(ft.name{ft_ix}) = model(:,ft_ix);
end

% Create Model Formula
reg_formula = strjoin(ft.name,' + ');
formula = ['ERP ~ ' reg_formula];

% Run Model: Single time point (mean in window) over channels
glm = cell(size(ch_list));
pvals = nan([numel(ft.name) numel(ch_list)]);
for ch_ix = 1:numel(ch_list)
    tbl.ERP = data(:,ch_ix);
    glm{ch_ix} = fitlme(tbl,formula);
    pvals(:,ch_ix) = glm{ch_ix}.Coefficients.pValue(2:end); % Skip intercept
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

%% Plot Regression results
fig_name = [SBJ_id '_' stat_id '_' an_id '_amp_fits'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
[n_rowcol,~] = fn_num_subplots(numel(ft.name));

for ft_ix = 1:numel(ft.name)
    subplot(n_rowcol(1),n_rowcol(2),ft_ix); hold on;
    
    % Compute correlation
    [r,p] = corrcoef(model(:,ft_ix),data);
    r = r(1,2); p = p(1,2);
    
    % Plot features
    scatter(model(:,ft_ix),data, 'o', 'k');
    
    % Plot linear fit
    coeff = polyfit(model(:,ft_ix),data,1);
    xbounds = get(gca,'XLim');
    xfudge = (xbounds(2)-xbounds(1))*0.1;
    xdat = [xbounds(1)+xfudge xbounds(2)-xfudge];
    ydat = coeff(1)*xdat + coeff(2);
    simple_fit = line(xdat,ydat);
    
    % Plot parameters
    legend(simple_fit,['corr r=' num2str(r,'%.3f') '; p=' num2str(p,'%.3f')],...
        'Location','best');
    title(['beta=' num2str(glm{1}.Coefficients.Estimate(ft_ix+1)) ...
        '; q=' num2str(qvals(ft_ix),'%.3f')]);
    xlabel([ft.name{ft_ix} '(' ft.chan{ft_ix} ', ' ft.cond{ft_ix} ') Amp (uV)']);
    ylabel('Mean Window Amp (uV)');
    set(gca,'FontSize',16);
end

% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/OB_feat/'];
    if ~exist(fig_dir,'dir') && save_fig
        mkdir(fig_dir);
    end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Save Results
stat_out_fname = [root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','glm','qvals','SBJs','time_vec','ch_list','reg_pk_time');

end
