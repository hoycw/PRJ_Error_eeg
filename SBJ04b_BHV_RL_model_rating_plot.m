function SBJ04b_BHV_RL_model_rating_plot(SBJ,proc_id,stat_id,varargin)
%% Plot single SBJ behavior with RL model fit and ratings
%   Scatter of single trial tolerance and outcomes
%   Scatter of block accuracy and average tolerance
%   Sigmoid from logistic regression fit on top
%   Scatter of single trial tolerance and ratings
%   Scatter of block averaged tolerance and ratings
%   Histograms (w/ t-test stats) for win vs. loss within easy/hard cond.
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       save_fig [0/1] - binary flag to save figure; default = 1
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figures

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
        elseif strcmp(varargin{v},'save_fig') && ischar(varargin{v+1})
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');     fig_vis = 'on'; end
if ~exist('fig_ftype','var');   fig_ftype = 'png'; end
if ~exist('save_fig','var');    save_fig = 1; end

%% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Plotting parameters
trl_sz      = 25;       % scatter size for single trial
mean_sz     = 100;      % scatter size for blocks
trl_rat_sz  = 15;       % scatter size for single trial ratings
mean_rat_sz = 75;      % scatter size for block ratings
sig_step    = 0.001;    % tolerance step size for plotting model fit

% Get model and condition parameters
[reg_lab, ~, ~, ~] = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.model_cond);

%% Load and Select Behavior
% Load data (should include betas from logistic regression)
load([SBJ_vars.dirs.proc SBJ '_model_' st.model_id '.mat']);

load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'],'bhv');
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

% Select Trails based on Conditions of Interest
full_cond_idx = fn_condition_index(cond_lab, bhv);
bhv_fields = fieldnames(bhv);
orig_n_trials = numel(bhv.trl_n);
for f_ix = 1:numel(bhv_fields)
    if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
        bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
    end
end
n_trials = numel(bhv.trl_n);
fprintf('\t%s: Loaded %d trials, kept %d for modeling...\n',SBJ,orig_n_trials,n_trials);

%% Compute Mean Accuracy, Ratings, and Tolerance per Block
% Adjust block numbers for EEG12
if strcmp(SBJ,'EEG12')
    blk5_starts = find(bhv.blk==5 & bhv.blk_trl_n==1);
    for trl_ix = 1:blk5_starts(2)-1
        bhv.blk(trl_ix) = bhv.blk(trl_ix)-4;
    end
end

% Compute block level accuracy, ratings, and tolerance
blk_ids = unique(bhv.blk);
blk_ez  = false(size(blk_ids));
blk_tol = nan(size(blk_ids));
blk_acc = nan(size(blk_ids));
blk_rat = nan(size(blk_ids));
for b_ix = 1:numel(blk_ids)
    blk_tol(b_ix) = mean(bhv.tol(bhv.blk==b_ix));
    blk_acc(b_ix) = mean(bhv.hit(bhv.blk==b_ix));
    blk_rat(b_ix) = nanmean(bhv.rating(bhv.blk==b_ix));
    if strcmp(unique(bhv.cond(bhv.blk==b_ix)),'easy')
        blk_ez(b_ix) = true;
    end
end

%% Compute Correlation between Subjective Ratings and Model Expected Value
% Get index for expected value and conditions
ev_idx     = strcmp(reg_lab,'EV');
ez_trl_idx = strcmp(bhv.cond,'easy');

ez_corr    = corrcoef(model(ez_trl_idx,ev_idx),bhv.rating(ez_trl_idx),'Rows','complete');
ez_corr    = ez_corr(1,2);
hd_corr    = corrcoef(model(~ez_trl_idx,ev_idx),bhv.rating(~ez_trl_idx),'Rows','complete');
hd_corr    = hd_corr(1,2);
total_corr = corrcoef(model(:,ev_idx),bhv.rating,'Rows','complete');
total_corr = total_corr(1,2);

%% Plot Tolerance vs. Outcome with Model Overlay
fig_name = [SBJ '_BHV_acc_' st.model_id '_pWin'];
figure('Name',fig_name,'Visible',fig_vis);
hold on;

% Plot Trial and Block Behavior: Easy
ez_trl = scatter(bhv.tol(ez_trl_idx), bhv.hit(ez_trl_idx), trl_sz,'k','filled');
ez_rat = scatter(bhv.tol(ez_trl_idx & bhv.hit==1), bhv.rating(ez_trl_idx & bhv.hit==1), trl_rat_sz,'g','*');
ez_rat = scatter(bhv.tol(ez_trl_idx & bhv.hit==0), bhv.rating(ez_trl_idx & bhv.hit==0), trl_rat_sz,'r','*');
ez_blk = scatter(blk_tol(blk_ez),blk_acc(blk_ez), mean_sz,'k','filled');
ez_rat_blk = scatter(blk_tol(blk_ez),blk_rat(blk_ez), mean_rat_sz,'b','filled','MarkerEdgeColor','k');

% Plot Trial and Block Behavior: Hard
hd_trl = scatter(bhv.tol(~ez_trl_idx), bhv.hit(~ez_trl_idx), trl_sz,'k','d');
hd_rat = scatter(bhv.tol(~ez_trl_idx & bhv.hit==1), bhv.rating(~ez_trl_idx & bhv.hit==1), trl_rat_sz,'g','+');
hd_rat = scatter(bhv.tol(~ez_trl_idx & bhv.hit==0), bhv.rating(~ez_trl_idx & bhv.hit==0), trl_rat_sz,'r','+');
hd_blk = scatter(blk_tol(~blk_ez),blk_acc(~blk_ez), mean_sz,'k','d');
hd_rat_blk = scatter(blk_tol(~blk_ez),blk_rat(~blk_ez), mean_rat_sz,'b','d','filled','MarkerEdgeColor','k');

% Reconstruct and plot Model Fit
sig_x = [0:sig_step:0.4];
sig_y = betas(1) + (sig_x * betas(2));
sig_y = 1 ./ (1+exp(-sig_y));
fit_line = line(sig_x,sig_y,'Color','k');

% Figure Parameters
xlabel('Tolerance (s)');
ylabel('Accuracy');
set(gca,'YLim',[0 1]);
title([SBJ ' Rating Correlations: r=' num2str(total_corr,'%.2f') ' (Easy=' ...
    num2str(ez_corr,'%.2f') '; Hard=' num2str(hd_corr,'%.2f') ')']);
ez_leg = ['Easy (n=' num2str(sum(ez_trl_idx)) '; mean = ' ...
    num2str(mean(bhv.hit(ez_trl_idx))*100,'%.1f') '%)'];
ez_rat_leg = ['Easy Ratings (mean = ' ...
    num2str(nanmean(bhv.rating(ez_trl_idx))*100,'%.1f') '%)'];
hd_leg = ['Hard (n=' num2str(sum(~ez_trl_idx)) '; mean = ' ...
    num2str(mean(bhv.hit(~ez_trl_idx))*100,'%.1f') '%)'];
hd_rat_leg = ['Hard Ratings (mean = ' ...
    num2str(nanmean(bhv.rating(~ez_trl_idx))*100,'%.1f') '%)'];
legend([ez_trl, ez_rat, hd_trl, hd_rat, fit_line],{ez_leg,ez_rat_leg,hd_leg,hd_rat_leg,'Model Fit'},'Location','southeast');
set(gca,'FontSize',14);

%% Save Model Scatter Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/rating_model_fits/' st.model_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    % Ensure vector graphics if saving
    if any(strcmp(fig_ftype,{'svg','eps'}))
        set(gcf, 'Renderer', 'painters');
    end
    saveas(gcf,fig_fname);
end

%% Statistics for win vs. loss ratings
% Statistics for win vs. loss within condition
[~, ez_pval] = ttest2(bhv.rating(ez_trl_idx & bhv.hit==0), bhv.rating(ez_trl_idx & bhv.hit==1));
[~, hd_pval] = ttest2(bhv.rating(~ez_trl_idx & bhv.hit==0), bhv.rating(~ez_trl_idx & bhv.hit==1));

% Plotting statistics
mean_ez_ls_rat = nanmean(bhv.rating(ez_trl_idx & bhv.hit==0));
mean_ez_wn_rat = nanmean(bhv.rating(ez_trl_idx & bhv.hit==1));
mean_hd_ls_rat = nanmean(bhv.rating(~ez_trl_idx & bhv.hit==0));
mean_hd_wn_rat = nanmean(bhv.rating(~ez_trl_idx & bhv.hit==1));

% Normalize ratings within condition
mean_ez_rat = nanmean(bhv.rating(ez_trl_idx));
std_ez_rat  = nanstd(bhv.rating(ez_trl_idx));
mean_hd_rat = nanmean(bhv.rating(~ez_trl_idx));
std_hd_rat  = nanstd(bhv.rating(~ez_trl_idx));

ez_wn_rat_norm = (bhv.rating(ez_trl_idx & bhv.hit==1) - mean_ez_rat)./std_ez_rat;
ez_ls_rat_norm = (bhv.rating(ez_trl_idx & bhv.hit==0) - mean_ez_rat)./std_ez_rat;
hd_wn_rat_norm = (bhv.rating(~ez_trl_idx & bhv.hit==1) - mean_hd_rat)./std_hd_rat;
hd_ls_rat_norm = (bhv.rating(~ez_trl_idx & bhv.hit==0) - mean_hd_rat)./std_hd_rat;

mean_ls_rat_norm = nanmean([ez_ls_rat_norm; hd_ls_rat_norm]);
mean_wn_rat_norm = nanmean([ez_wn_rat_norm; hd_wn_rat_norm]);

% Stats: win vs. loss combined after normalizing for mean within easy/hard
[~, norm_pval] = ttest2([ez_wn_rat_norm; hd_wn_rat_norm], [ez_ls_rat_norm; hd_ls_rat_norm]);

%% Plot rating histograms by condition and outcome
% Set Up histogram plots
n_bins = 20;
fig_name = [SBJ '_BHV_ratings_WinVsLoss_cond'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1.0 0.5],'Visible',fig_vis);

% Plot Easy histograms
subplot(1,3,1); hold on;
histogram(bhv.rating(ez_trl_idx & bhv.hit==1),n_bins,'FaceColor','g','FaceAlpha',0.3);
histogram(bhv.rating(ez_trl_idx & bhv.hit==0),n_bins,'FaceColor','r','FaceAlpha',0.3);
ez_ls_line = line([mean_ez_ls_rat mean_ez_ls_rat], ylim, 'Color', 'r', 'LineWidth', 3);
ez_wn_line = line([mean_ez_wn_rat mean_ez_wn_rat], ylim, 'Color', 'g', 'LineWidth', 3);
xlabel('Subjective Win % Rating');
xlim([0 1]);
ylabel('# Trials');
legend([ez_wn_line ez_ls_line],{['Easy Wins: n=' num2str(sum(~isnan(bhv.rating(ez_trl_idx & bhv.hit==1)))) ...
    '; mean=' num2str(mean_ez_wn_rat,'%.3f')], ['Easy Losses: n=' num2str(sum(~isnan(bhv.rating(ez_trl_idx & bhv.hit==0))))...
    '; mean=' num2str(mean_ez_ls_rat,'%.3f')]},'Location','northwest');
title(['Easy Win vs. Loss: p = ', num2str(ez_pval,'%.3f')]);
set(gca,'FontSize',14);

% Plot Hard histograms
subplot(1,3,2); hold on;
histogram(bhv.rating(~ez_trl_idx & bhv.hit==0),n_bins,'FaceColor','r','FaceAlpha',0.3);
histogram(bhv.rating(~ez_trl_idx & bhv.hit==1),n_bins,'FaceColor','g','FaceAlpha',0.3);
hd_ls_line = line([mean_hd_ls_rat mean_hd_ls_rat], ylim, 'Color', 'r', 'LineWidth', 3);
hd_wn_line = line([mean_hd_wn_rat mean_hd_wn_rat], ylim, 'Color', 'g', 'LineWidth', 3);
xlabel('Subjective Win % Rating');
xlim([0 1]);
ylabel('# Trials');
legend([hd_wn_line hd_ls_line],{['Hard Wins: n=' num2str(sum(~isnan(bhv.rating(~ez_trl_idx & bhv.hit==1)))) ...
    '; mean=' num2str(mean_hd_wn_rat,'%.3f')], ['Hard Losses: n=' num2str(sum(~isnan(bhv.rating(~ez_trl_idx & bhv.hit==0))))...
    '; mean=' num2str(mean_hd_ls_rat,'%.3f')]},'Location','northeast');
title(['Hard Win vs. Loss: p = ', num2str(hd_pval,'%.3f')]);
set(gca,'FontSize',14);

% Plot Normalized histograms
subplot(1,3,3); hold on;
histogram([ez_ls_rat_norm; hd_ls_rat_norm],n_bins,'FaceColor','r','FaceAlpha',0.3);
histogram([ez_wn_rat_norm; hd_wn_rat_norm],n_bins,'FaceColor','g','FaceAlpha',0.3);
norm_ls_line = line([mean_ls_rat_norm mean_ls_rat_norm], ylim, 'Color', 'r', 'LineWidth', 3);
norm_wn_line = line([mean_wn_rat_norm mean_wn_rat_norm], ylim, 'Color', 'g', 'LineWidth', 3);
xlabel('Z-Scored Subjective Win % Rating');
%xlim([0 1]);
ylabel('# Trials');
legend([norm_wn_line norm_ls_line],{['Wins: n=' num2str(sum(~isnan(bhv.rating(bhv.hit==1)))) ...
    '; mean=' num2str(mean_wn_rat_norm,'%.3f')], ['Losses: n=' num2str(sum(~isnan(bhv.rating(bhv.hit==0))))...
    '; mean=' num2str(mean_ls_rat_norm,'%.3f')]},'Location','best');
title(['Normalized Win vs. Loss: p = ', num2str(norm_pval,'%.3f')]);
set(gca,'FontSize',14);

%% Save Ratings Histogram Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/rating_WvsL/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    % Ensure vector graphics if saving
    if any(strcmp(fig_ftype,{'svg','eps'}))
        set(gcf, 'Renderer', 'painters');
    end
    saveas(gcf,fig_fname);
end

end
