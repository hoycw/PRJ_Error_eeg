function SBJ04b_BHV_RL_model_plot(SBJ,proc_id,stat_id,varargin)
% Plot single SBJ tolerance, accuracy, outcomes with RL model fit
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
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
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('save_fig','var'); save_fig = 1; end

%% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, cond_colors, ~, ~] = fn_condition_label_styles(st.trial_cond{1});
[eh_lab, ~, eh_colors, ~, ~] = fn_condition_label_styles('Dif');

%% Load and Select Behavior
% Load data
load([SBJ_vars.dirs.proc SBJ '_model_' model_id '.mat']);

load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'],'bhv');
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

% Select Conditions of Interest
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

%% Compute Mean Accuracy and Tolerance per Block
% Adjust block numbers for EEG12
if strcmp(SBJ,'EEG12')
    blk5_starts = find(bhv.blk==5 & bhv.blk_trl_n==1);
    for trl_ix = 1:blk5_starts(2)-1
        bhv.blk(trl_ix) = bhv.blk(trl_ix)-4;
    end
end

% Compute stats
blk_ids = unique(bhv.blk);
blk_ez  = false(size(blk_ids));
blk_tol = nan(size(blk_ids));
blk_acc = nan(size(blk_ids));
for b_ix = 1:numel(blk_ids)
    blk_tol(b_ix) = mean(bhv.tol(bhv.blk==b_ix));
    blk_acc(b_ix) = mean(bhv.hit(bhv.blk==b_ix));
    if strcmp(unique(bhv.cond(bhv.blk==b_ix)),'easy')
        blk_ez(b_ix) = true;
    end
end

%% Plot Tolerance vs. Outcome
fig_name = [SBJ '_BHV_acc_' model_id '_pWin'];
figure('Name',fig_name,'Visible',fig_vis);
hold on;
trl_sz = 25;
mean_sz = 100;
sig_step = 0.001;

% Plot Trial and Block Behavior: Easy
ez_trl_idx = strcmp(bhv.cond,'easy');
ez_trl = scatter(bhv.tol(ez_trl_idx), bhv.hit(ez_trl_idx), trl_sz,'k','filled');
ez_blk = scatter(blk_tol(blk_ez),blk_acc(blk_ez), mean_sz,'k','filled');

% Plot Trial and Block Behavior: Hard
hd_trl = scatter(bhv.tol(~ez_trl_idx), bhv.hit(~ez_trl_idx), trl_sz,'k','d');
hd_blk = scatter(blk_tol(~blk_ez),blk_acc(~blk_ez), mean_sz,'k','d');

% Plot Model Fit
sig_x = [0:sig_step:0.4];
sig_y = betas(1) + (sig_x * betas(2));
sig_y = 1 ./ (1+exp(-sig_y));
fit_line = line(sig_x,sig_y,'Color','k');

% Figure Parameters
xlabel('Tolerance (s)');
ylabel('Accuracy');
set(gca,'YLim',[0 1]);
title(SBJ);
ez_leg = ['Easy (n=' num2str(sum(ez_trl_idx)) '; mean = ' ...
    num2str(mean(bhv.hit(ez_trl_idx))*100,'%.1f') '%)'];
hd_leg = ['Hard (n=' num2str(sum(~ez_trl_idx)) '; mean = ' ...
    num2str(mean(bhv.hit(~ez_trl_idx))*100,'%.1f') '%)'];
legend([ez_trl, hd_trl, fit_line],{ez_leg,hd_leg,'Model Fit'},'Location','southeast');
set(gca,'FontSize',14);

%% Save Figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/model_fits/' model_id '/'];
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
