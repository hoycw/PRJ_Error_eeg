function BHV05_grp_rating_stats(SBJ_id,proc_id,stat_id,varargin)
%% Group stats for subjective ratings
%   (1) correlation of rating vs. pWin regressor
%   (2) learning: correlation values across blocks
%   (3) t-test for ratings greater for wins than losses
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       save_fig [0/1] - binary flag to save figure; default = 1
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   (1) Scatter plot of rating vs. pWin for easy, hard, all trials
%   (2) Line plot of correlation block-by-block for easy, hard, all trials
%   (3) Histograms of z-scored Ratings for win vs. loss (easy, hard, all)

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
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[reg_lab, ~, ~, ~] = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.model_cond);
ev_ix = strcmp(reg_lab,'EV');

win_color  = [51 160 44]./256;
loss_color = [227 26 28]./256;

%% Load Behavior
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
rate_rt_thresh= zeros([numel(SBJs) 1]);
for s = 1:numel(SBJs)
    % Load data
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    bhvs{s} = tmp.bhv;
    
    % Select Conditions of Interest
    full_cond_idx{s} = fn_condition_index(cond_lab, bhvs{s});
    
    % Select only rating trials
    non_rating_idx = isnan(bhvs{s}.rating);
    
    % Remove bad ratings (RT outliers, time outs)
    time_out_idx = bhvs{s}.rating_time_out==1;
    rate_rt_thresh(s) = nanmean(bhvs{s}.rating_rt(~time_out_idx))+nanstd(bhvs{s}.rating_rt(~time_out_idx))*3;
    rt_out_idx = bhvs{s}.rating_rt >= rate_rt_thresh(s);
    full_cond_idx{s}(non_rating_idx | time_out_idx | rt_out_idx) = 0;
    fprintf('\t%s removed %i trials: %i non-ratings, %i rt outliers (mean=%.2f), %i time outs\n',SBJs{s},...
        sum(non_rating_idx)+sum(rt_out_idx)+sum(time_out_idx),sum(non_rating_idx),sum(rt_out_idx),...
        mean(bhvs{s}.rating_rt(rt_out_idx)),sum(time_out_idx));
    
    % Plot Rating RT dist
%     histogram(bhvs{s}.rating_rt); line([rate_rt_thresh(s) rate_rt_thresh(s)],ylim); pause;
    
    % Remove trials
    bhv_fields = fieldnames(bhvs{s});
    orig_n_trials = numel(bhvs{s}.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhvs{s}.(bhv_fields{f_ix}))==orig_n_trials
            bhvs{s}.(bhv_fields{f_ix}) = bhvs{s}.(bhv_fields{f_ix})(full_cond_idx{s}~=0);
        end
    end
    n_trials(s) = numel(bhvs{s}.trl_n);
    
    clear tmp time_out_idx rt_thresh rt_out_idx
end

%% Load Data and Build Model
data = struct;
data.pWin    = zeros([sum(n_trials) 1]);
data.rating  = zeros([sum(n_trials) 1]);
data.zrating = zeros([sum(n_trials) 1]);
data.sbj     = zeros([sum(n_trials) 1]);
data.blk     = zeros([sum(n_trials) 1]);
data.blk_acc = zeros([sum(n_trials) 1]);
data.blk_rat = zeros([sum(n_trials) 1]);
data.ez      = zeros([sum(n_trials) 1]);
data.hit     = zeros([sum(n_trials) 1]);
for s = 1:numel(SBJs)
    % Track SBJ in design matrix
    if s==1
        data.sbj(1:n_trials(s)) = s*ones([n_trials(s) 1]);
    else
        data.sbj(sum(n_trials(1:s-1))+1:sum(n_trials(1:s))) = s*ones([n_trials(s) 1]);
    end
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' st.model_id '.mat']);
    
    % Z-score rating within condition
    sbj_ez_idx = strcmp(bhvs{s}.cond,'easy');
    sbj_zrate = nan(size(sbj_ez_idx));
    zrate_ez = (bhvs{s}.rating(sbj_ez_idx)-nanmean(bhvs{s}.rating(sbj_ez_idx)))./nanstd(bhvs{s}.rating(sbj_ez_idx));
    zrate_hd = (bhvs{s}.rating(~sbj_ez_idx)-nanmean(bhvs{s}.rating(~sbj_ez_idx)))./nanstd(bhvs{s}.rating(~sbj_ez_idx));
    sbj_zrate(sbj_ez_idx) = zrate_ez;
    sbj_zrate(~sbj_ez_idx) = zrate_hd;
    
    % Compute block accuracy
    blk_ids = unique(bhvs{s}.blk);
    blk_acc = nan(size(bhvs{s}.blk));
    blk_rat = nan(size(bhvs{s}.blk));
    for b_ix = 1:numel(blk_ids)
        tmp_acc = mean(bhvs{s}.hit(bhvs{s}.blk==b_ix));
        tmp_rat = nanmean(bhvs{s}.rating(bhvs{s}.blk==b_ix));
        blk_acc(bhvs{s}.blk==b_ix) = tmp_acc*ones([sum(bhvs{s}.blk==b_ix) 1]);
        blk_rat(bhvs{s}.blk==b_ix) = tmp_rat*ones([sum(bhvs{s}.blk==b_ix) 1]);
    end
    
    % Compile data
    data.pWin(data.sbj==s)    = (tmp.model(full_cond_idx{s}~=0,ev_ix) + 1)/2;   % Convert Expected Value back to Win Probability
    data.rating(data.sbj==s)  = bhvs{s}.rating;
    data.zrating(data.sbj==s) = sbj_zrate;
    data.blk(data.sbj==s)     = bhvs{s}.blk;
    data.blk_acc(data.sbj==s) = blk_acc;
    data.blk_rat(data.sbj==s) = blk_rat;
    data.ez(data.sbj==s)      = sbj_ez_idx;
    data.hit(data.sbj==s)     = bhvs{s}.hit;
    
    clear tmp sbj_zrate sbj_ez_idx zrate_ez zrate_hd
end

%% Compute Correlation between Subjective Ratings and Model Expected Value
[tmp_corr,tmp_p] = corrcoef(data.pWin(data.ez==1),data.rating(data.ez==1),'Rows','complete');
ez_corr = tmp_corr(1,2); ez_pval = tmp_p(1,2);

[tmp_corr,tmp_p] = corrcoef(data.pWin(data.ez==0),data.rating(data.ez==0),'Rows','complete');
hd_corr = tmp_corr(1,2); hd_pval = tmp_p(1,2);

[tmp_corr,tmp_p] = corrcoef(data.pWin,data.rating,'Rows','complete');
full_corr = tmp_corr(1,2); full_pval = tmp_p(1,2);

[tmp_corr,tmp_p] = corrcoef(data.blk_acc,data.rating,'Rows','complete');
acc_corr = tmp_corr(1,2); acc_pval = tmp_p(1,2);

%% Plot pWin and Accuracy Scatters
fig_name = [SBJ_id '_BHV_ratings_scatter_corr'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1.0 1.0],'Visible',fig_vis);

trl_rat_sz = 15;
pWin_fudge = 0.05;
sig_step   = 0.001;    % tolerance step size for plotting model fit
lin_fit_x  = [0.1:sig_step:0.9];

% Plot pWin: Easy Trials
subplot(2,2,1); hold on;
win_scat  = scatter(data.rating(data.ez==1 & data.hit==1),data.pWin(data.ez==1 & data.hit==1),trl_rat_sz,win_color,'o');
loss_scat = scatter(data.rating(data.ez==1 & data.hit==0),data.pWin(data.ez==1 & data.hit==0),trl_rat_sz,loss_color,'*');

lin_fit = polyfit(data.rating(data.ez==1),data.pWin(data.ez==1),1);
lin_fit_y = lin_fit(1)*lin_fit_x + lin_fit(2);
plot(lin_fit_x,lin_fit_y,'k','LineWidth',2);

xlabel('Subjective Rating');
ylabel('Model Win %');
xlim([0 1]);
ylim([0 1]);
legend([win_scat loss_scat],{'Win Ratings','Loss Ratings'},'Location','southwest');
title(['Easy Rating Correlations: r=' num2str(ez_corr,'%.2f') '; p=' num2str(ez_pval,'%.3f')]);
set(gca,'FontSize',14);

% Plot pWin: Hard Trials
subplot(2,2,2); hold on;
win_scat  = scatter(data.rating(data.ez==0 & data.hit==1),data.pWin(data.ez==0 & data.hit==1),trl_rat_sz,win_color,'o');
loss_scat = scatter(data.rating(data.ez==0 & data.hit==0),data.pWin(data.ez==0 & data.hit==0),trl_rat_sz,loss_color,'*');

lin_fit = polyfit(data.rating(data.ez==0),data.pWin(data.ez==0),1);
lin_fit_y = lin_fit(1)*lin_fit_x + lin_fit(2);
plot(lin_fit_x,lin_fit_y,'k','LineWidth',2);

xlabel('Subjective Rating');
ylabel('Model Win %');
xlim([0 1]);
ylim([0 1]);
legend([win_scat loss_scat],{'Win Ratings','Loss Ratings'},'Location','northeast');
title(['Hard Rating Correlations: r=' num2str(hd_corr,'%.2f') '; p=' num2str(hd_pval,'%.3f')]);
set(gca,'FontSize',14);

% Plot pWin: All Trials
subplot(2,2,3); hold on;
scatter(data.rating(data.hit==1),data.pWin(data.hit==1),trl_rat_sz,win_color,'o');
scatter(data.rating(data.hit==0),data.pWin(data.hit==0),trl_rat_sz,loss_color,'*');

lin_fit = polyfit(data.rating,data.pWin,1);
lin_fit_y = lin_fit(1)*lin_fit_x + lin_fit(2);
plot(lin_fit_x,lin_fit_y,'k','LineWidth',2);

xlabel('Subjective Rating');
ylabel('Model Win %');
xlim([0 1]);
ylim([0 1]);
title(['All Trials Rating Correlations: r=' num2str(full_corr,'%.2f') '; p=' num2str(full_pval,'%.3f')]);
set(gca,'FontSize',14);

% Plot Block Accuracy: All Trials
subplot(2,2,4); hold on;
scatter(data.rating(data.hit==1),data.blk_acc(data.hit==1),trl_rat_sz,win_color,'o');
scatter(data.rating(data.hit==0),data.blk_acc(data.hit==0),trl_rat_sz,loss_color,'*');

lin_fit = polyfit(data.rating,data.blk_acc,1);
lin_fit_y = lin_fit(1)*lin_fit_x + lin_fit(2);
plot(lin_fit_x,lin_fit_y,'k');

xlabel('Subjective Rating');
ylabel('Block Accuracy');
xlim([0 1]);
ylim([0 1]);
title(['All Trials Correlations: r=' num2str(acc_corr,'%.2f') '; p=' num2str(acc_pval,'%.3f')]);
set(gca,'FontSize',14);

% Save Correlation Scatter Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/rating_pWin_scatter/' st.model_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compute Block-by-block correlations
blk_ids = unique(data.blk);
blk_corrs = nan([numel(blk_ids) 3]);
blk_pvals = nan([numel(blk_ids) 3]);
for b_ix = 1:numel(blk_ids)
    % Easy
    [tmp_corr,tmp_p] = corrcoef(data.pWin(data.blk==b_ix & data.ez==1),data.rating(data.blk==b_ix & data.ez==1),'Rows','complete');
    blk_corrs(b_ix,1) = tmp_corr(1,2); blk_pvals(b_ix,1) = tmp_p(1,2);
    
    % Hard
    [tmp_corr,tmp_p] = corrcoef(data.pWin(data.blk==b_ix & data.ez==0),data.rating(data.blk==b_ix & data.ez==0),'Rows','complete');
    blk_corrs(b_ix,2) = tmp_corr(1,2); blk_pvals(b_ix,2) = tmp_p(1,2);
    
    % All trials
    [tmp_corr,tmp_p] = corrcoef(data.pWin(data.blk==b_ix),data.rating(data.blk==b_ix),'Rows','complete');
    blk_corrs(b_ix,3) = tmp_corr(1,2); blk_pvals(b_ix,3) = tmp_p(1,2);
end

%% Plot Block-by-Block Correlations
fig_name = [SBJ_id '_BHV_ratings_block_corr'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1.0 0.5],'Visible',fig_vis);

[eh_lab,eh_names,eh_colors,~,~] = fn_condition_label_styles('Dif');
hold on;

ez_line = line(1:numel(blk_ids),blk_corrs(:,1),'Color',eh_colors{1});
hd_line = line(1:numel(blk_ids),blk_corrs(:,2),'Color',eh_colors{2});
full_line = line(1:numel(blk_ids),blk_corrs(:,3),'Color','k');

xlim([0 numel(blk_ids)+1]);
xticks([1:numel(blk_ids)]);
xlabel('Block');
ylabel('Correlation (r)');
legend([ez_line hd_line full_line],{'Easy','Hard','All Trials'});
set(gca,'FontSize',14);

% Save Correlation Scatter Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/rating_block_corr/' st.model_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot and Print Bias between pWin and Ratings
% Average pWin and Ratings within easy and hard
blk_pwin = nan([numel(SBJs) 2]);
blk_rate = nan([numel(SBJs) 2]);
for s = 1:numel(SBJs)
    for ez_ix = [0 1]
        blk_pwin(s,ez_ix+1) = mean(data.pWin(data.sbj==s & data.ez==ez_ix));
        blk_rate(s,ez_ix+1) = mean(data.rating(data.sbj==s & data.ez==ez_ix));
    end
end
pwin_rate_diff = blk_pwin-blk_rate;

fig_name = [SBJ_id '_BHV_ratings_pWin_bias'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);

violins = violinplot(pwin_rate_diff, {'Hard','Easy'}, 'ShowData', true, 'ShowMean', true, 'ViolinAlpha', 0.3);
violins(1).ViolinColor = eh_colors{2};
violins(2).ViolinColor = eh_colors{1};
ylabel('Bias: (pWin - Rating)');
legend([violins(1).ViolinPlot violins(2).ViolinPlot],...
    {['Hard mean = ' num2str(mean(pwin_rate_diff(:,1)),'%.2f')],...
    ['Easy mean = ' num2str(mean(pwin_rate_diff(:,2)),'%.2f')]},'Location','best');
set(gca,'FontSize',14);

% Save Ratings Histogram Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/rating_pWin_bias/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Statistics for win vs. loss ratings
% Statistics for win vs. loss within condition
[~, ez_pval] = ttest2(data.zrating(data.ez==1 & data.hit==0), data.zrating(data.ez==1 & data.hit==1));
[~, hd_pval] = ttest2(data.zrating(data.ez==0 & data.hit==0), data.zrating(data.ez==0 & data.hit==1));
[~, full_pval] = ttest2(data.zrating(data.hit==0), data.zrating(data.hit==1));

% Plotting statistics
mean_ez_ls_rat = nanmean(data.zrating(data.ez==1 & data.hit==0));
mean_ez_wn_rat = nanmean(data.zrating(data.ez==1 & data.hit==1));
mean_hd_ls_rat = nanmean(data.zrating(data.ez==0 & data.hit==0));
mean_hd_wn_rat = nanmean(data.zrating(data.ez==0 & data.hit==1));
mean_fl_ls_rat = nanmean(data.zrating(data.hit==0));
mean_fl_wn_rat = nanmean(data.zrating(data.hit==1));

%% Plot rating histograms by condition and outcome
% Set Up histogram plots
n_bins = 20;
fig_name = [SBJ_id '_BHV_ratings_WinVsLoss_cond'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1.0 0.5],'Visible',fig_vis);

% Plot Easy histograms
subplot(1,3,1); hold on;
histogram(data.zrating(data.ez==1 & data.hit==1),n_bins,'FaceColor',win_color,'FaceAlpha',0.3);
histogram(data.zrating(data.ez==1 & data.hit==0),n_bins,'FaceColor',loss_color,'FaceAlpha',0.3);
ez_ls_line = line([mean_ez_ls_rat mean_ez_ls_rat], ylim, 'Color', loss_color, 'LineWidth', 3);
ez_wn_line = line([mean_ez_wn_rat mean_ez_wn_rat], ylim, 'Color', win_color, 'LineWidth', 3);
xlabel('Z-Scored Subjective Win % Rating');
ylabel('# Trials');
legend([ez_wn_line ez_ls_line],{['Easy Wins: n=' num2str(sum(data.ez==1 & data.hit==1)) ...
    '; mean=' num2str(mean_ez_wn_rat,'%.3f')], ['Easy Losses: n=' num2str(sum(data.ez==1 & data.hit==0))...
    '; mean=' num2str(mean_ez_ls_rat,'%.3f')]},'Location','northwest');
title(['Easy Win vs. Loss: p = ', num2str(ez_pval,'%.3f')]);
set(gca,'FontSize',14);

% Plot Hard histograms
subplot(1,3,2); hold on;
histogram(data.zrating(data.ez==0 & data.hit==1),n_bins,'FaceColor',win_color,'FaceAlpha',0.3);
histogram(data.zrating(data.ez==0 & data.hit==0),n_bins,'FaceColor',loss_color,'FaceAlpha',0.3);
hd_ls_line = line([mean_hd_ls_rat mean_hd_ls_rat], ylim, 'Color', loss_color, 'LineWidth', 3);
hd_wn_line = line([mean_hd_wn_rat mean_hd_wn_rat], ylim, 'Color', win_color, 'LineWidth', 3);
xlabel('Z-Scored Subjective Win % Rating');
ylabel('# Trials');
legend([hd_wn_line hd_ls_line],{['Hard Wins: n=' num2str(sum(data.ez==0 & data.hit==1)) ...
    '; mean=' num2str(mean_hd_wn_rat,'%.3f')], ['Hard Losses: n=' num2str(sum(data.ez==0 & data.hit==0))...
    '; mean=' num2str(mean_hd_ls_rat,'%.3f')]},'Location','northwest');
title(['Hard Win vs. Loss: p = ', num2str(hd_pval,'%.3f')]);
set(gca,'FontSize',14);

% Plot All Trial histograms
subplot(1,3,3); hold on;
histogram(data.zrating(data.hit==1),n_bins,'FaceColor',win_color,'FaceAlpha',0.3);
histogram(data.zrating(data.hit==0),n_bins,'FaceColor',loss_color,'FaceAlpha',0.3);
norm_ls_line = line([mean_fl_ls_rat mean_fl_ls_rat], ylim, 'Color', loss_color, 'LineWidth', 3);
norm_wn_line = line([mean_fl_wn_rat mean_fl_wn_rat], ylim, 'Color', win_color, 'LineWidth', 3);
xlabel('Z-Scored Subjective Win % Rating');
%xlim([0 1]);
ylabel('# Trials');
legend([norm_wn_line norm_ls_line],{['Wins: n=' num2str(sum(data.hit==1)) ...
    '; mean=' num2str(mean_fl_wn_rat,'%.3f')], ['Losses: n=' num2str(sum(data.hit==0))...
    '; mean=' num2str(mean_fl_ls_rat,'%.3f')]},'Location','best');
title(['Normalized Win vs. Loss: p = ', num2str(full_pval,'%.3f')]);
set(gca,'FontSize',14);

% Save Ratings Histogram Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/rating_WvsL/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
