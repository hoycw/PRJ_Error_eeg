function SBJ06d_OB_TT_ERP_grp_stats_corr_ts(SBJ_id,tt_proc_id,ob_proc_id,stat_id,plt_id,varargin)
%% Correlate OB ERP features with TT ERP over time per condition
%   Oddball ERPs only yield one feature per SBJ, so predict individual OB-TT
%   feature-condition pairs separately
% COMPUTATIONS:
%   Load OB and TT ERP features
%       Optional: z-score regressors across group
%   Correlate each OB feature with each TT condition at every time point
%   Correct for multiple comparisons (FDR for regressors)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   tt_proc_id [str] - ID of target time preprocessing pipeline
%   ob_proc_id [str] - ID of oddball preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%       st.model   = feat_id for OB ERP features
%       st.measure = 'ts'
%       st.an_id = TT ERPs to predict
% OUTPUTS:
%   corrs [cell array] -  [n_OB_feat, n_TT_cond, n_time] r values
%   pvals [float array] - [n_OB_feat, n_TT_cond, n_time] p values
%   qvals [float array] - [n_OB_feat, n_TT_cond, n_time] p values adjusted for multiple comparisons 
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)

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
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/OB_TT_feat/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.an_style,'corr'); error('This script is for correlation'); end
if ~strcmp(st.measure,'ts'); error('This script is for analysis across time!'); end

% OB Feature Parameters
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' st.model_lab '_vars.m'];
eval(feat_vars_cmd);
if ~any(strcmp(ft.grp_id,{'rare','Odd','Tar'})); error('Features should be oddball conditions!'); end

% Plot Parameters
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, cond_names, cond_colors, cond_styles, ~] = fn_condition_label_styles(st.stat_cond);

%% Load Oddball ERP Features
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'SBJs');
if numel(SBJs)~=numel(tmp.SBJs) || ~all(strcmp(SBJs,tmp.SBJs)); error('SBJ mismatch!'); end
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'ft_amp','ft_times');

% Z-score feature predictors
model = ft_amp;
if st.z_reg
    error('why z-score OB ERP features for correlation?');
    % model = zscore(model);
end

%% Load Behavior and Select Conditions
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

%% Load TT Data and Compute ERPs
cfgs  = [];
cfgs.latency = st.stat_lim;
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',st.an_id,'.mat'],'roi');
    if numel(roi.label)>1; error('Only run for one channel!'); end
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_roi = ft_selectdata(cfgs, roi);
    
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = st_roi.time{1};
        erps     = nan([numel(cond_lab) numel(SBJs) numel(time_vec)]);
    end
    
    % Compute ERPs within condition
    cond_idx = fn_condition_index(cond_lab, bhvs{s});
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        trials = nan([numel(cond_trial_ix) numel(time_vec)]);
        for trl_ix = 1:numel(cond_trial_ix)
            trials(trl_ix,:) = st_roi.trial{cond_trial_ix(trl_ix)};
        end
        erps(cond_ix,s,:) = mean(trials,1);
    end
    
    clear tmp roi ft_roi trials cond_idx cond_trial_ix
end

%% Run Correlation Over Time
fprintf('========================== Running Stats ==========================\n');
corrs = nan([numel(ft.name) numel(cond_lab) numel(time_vec)]);
pvals = nan([numel(ft.name) numel(cond_lab) numel(time_vec)]);
for ft_ix = 1:numel(ft.name)
    for cond_ix = 1:numel(cond_lab)
        for t_ix = 1:numel(time_vec)
            [r,p] = corrcoef(model(:,ft_ix),erps(cond_ix,:,t_ix),'Rows','complete');
            corrs(ft_ix,cond_ix,t_ix) = r(1,2);
            pvals(ft_ix,cond_ix,t_ix) = p(1,2);
        end
    end
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2)*size(pvals,3) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2) size(pvals,3)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

%% Plot Results
% Compute means and variance
plot_means = NaN([numel(cond_lab) numel(time_vec)]);
sems  = NaN([numel(cond_lab) numel(time_vec)]);
for cond_ix = 1:numel(cond_lab)
    plot_means(cond_ix,:) = squeeze(mean(erps(cond_ix,:,:),2));
    sems(cond_ix,:) = squeeze(std(erps(cond_ix,:,:),[],2))./sqrt(numel(SBJs))';
end

% Plot ERPs and condition correlations per OB feature
for ft_ix = 1:numel(ft.name)
    fig_name = [SBJ_id '_' stat_id '_' ft.name{ft_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.7 1],'Visible',fig_vis);
    
    % Find significant time periods
    sig_cond   = false(size(cond_lab));
    sig_chunks = cell(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        if any(pvals(ft_ix,cond_ix,:) <= st.alpha)
            sig_cond(cond_ix) = true;
            
            % Find consecutive chunks of (non-)significant coefficients
            sig_chunks{cond_ix} = fn_find_chunks(squeeze(pvals(ft_ix,cond_ix,:))<=st.alpha);
            % Remove non-significant chunks
            sig_chunks{cond_ix}(squeeze(pvals(ft_ix,cond_ix,sig_chunks{cond_ix}(:,1))>st.alpha),:) = [];
        end
    end
    
    %% Plot ERPs
    axes = gobjects([2 1]);
    subplot(2,1,1);
    axes(1) = gca; hold on;
    
    % Plot ERP Means (and variance)
    cond_lines = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(plt.evnt_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        cond_lines{cond_ix} = shadedErrorBar(time_vec, plot_means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(cond_ix) = cond_lines{cond_ix}.mainLine;
    end
    
    % Fix y limits to be consistent across channels
    ylims = [-15 30];
    %     if any(strcmp(SBJ_id,{'good1','goodall'}))
    %     else
    %         ylims = ylim;
    %     end
    
    % Plot Feedback onset
    main_lines(numel(cond_lab)+1) = line(...
        [0 0],[-15 30],...%ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    leg_lab = [cond_names plt.evnt_lab];
    
    % Axes and Labels
    axes(1).YLabel.String = 'uV';
    axes(1).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(1).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(1).XLabel.String = 'Time (s)';
    axes(1).Title.String  = ['TT ERPs (n = ' num2str(numel(SBJs)) ')'];
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    % Fix y limits to be consistent across channels
    if any(strcmp(SBJ_id,{'good1','goodall'}))
        ylims = [-15 30];
    else
        ylims = ylim;
    end
    set(gca,'FontSize',16);
    axes(1).YLim = ylims;
    
    %% Plot correlation values
    subplot(2,1,2);
    axes(2) = gca; hold on;
    
    corr_lines = gobjects([numel(cond_lab)+numel(plt.evnt_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        corr_lines(cond_ix) = line(time_vec, squeeze(corrs(ft_ix,cond_ix,:)),...
            'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix});
    
        % Plot Significance by bolding time series
        for sig_ix = 1:size(sig_chunks{cond_ix},1)
            % Assume: strcmp(plt.sig_type,'bold')
            sig_times = time_vec(sig_chunks{cond_ix}(sig_ix,1):sig_chunks{cond_ix}(sig_ix,2));
            line(sig_times,squeeze(corrs(ft_ix,cond_ix,sig_chunks{cond_ix}(sig_ix,1):sig_chunks{cond_ix}(sig_ix,2))),...
                'Color',cond_colors{cond_ix},'LineStyle',cond_styles{cond_ix},...
                'LineWidth',plt.sig_width);
        end
    end
    
    % Plot Feedback onset
    corr_lines(numel(cond_lab)+1) = line(...
        [0 0],ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
        'LineStyle',plt.evnt_styles{1});
    
    % Axes and Labels
    axes(2).YLim          = [-1 1];
    axes(2).YLabel.String = 'Correlation (r)';
    axes(2).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(2).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(2).XLabel.String = 'Time (s)';
    axes(2).Title.String  = [ft.name{ft_ix} ': OB-TT Correlations'];
    if plt.legend
        legend(corr_lines,[cond_names plt.evnt_lab],'Location',plt.legend_loc);
    end
    set(gca,'FontSize',16);
    
    %% Report peak stats per regressor
    % Prints largest model coefficient, time point, and q value
    for cond_ix = 1:numel(cond_lab)
        max_tmp = max(corrs(ft_ix,cond_ix,:));
        min_tmp = min(corrs(ft_ix,cond_ix,:));
        if abs(max_tmp) > abs(min_tmp)
            [max_corr, max_t_ix] = max(corrs(ft_ix,cond_ix,:));
        else
            [max_corr, max_t_ix] = min(corrs(ft_ix,cond_ix,:));
        end
        fprintf('%s: %s max corr = %.03f at %.03f; p = %.03f\n',ft.name{ft_ix},...
            cond_lab{cond_ix},max_corr,time_vec(max_t_ix),pvals(ft_ix,cond_ix,max_t_ix));
    end
    
    %% Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        % Ensure vector graphics if saving
        if any(strcmp(fig_ftype,{'svg','eps'}))
            set(gcf, 'Renderer', 'painters');
        end
        saveas(gcf,fig_fname);
    end
end

%% Save Results
stat_out_fname = [root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','corrs','pvals','qvals','SBJs');

end
