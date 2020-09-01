function SBJ04b_RL_model_plot_grp_predictions(SBJ_id,proc_id,stat_id,plt_id,save_fig,varargin)
%% Plot model predictors across group averaged within condition
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   proc_id [str] - ID of preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%   plt_id [str] - ID of the plotting parameters to use
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
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
if ischar(save_fig);            save_fig = str2num(save_fig); end

%% Load Data 
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, reg_names, reg_colors, reg_styles, reg_mrkrs] = fn_regressor_label_styles(st.model_lab);
[cond_lab, cond_names, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});
ez_idx = ~cellfun(@isempty,strfind(cond_lab,'Ez'));

%% Load Model and Comute Means
model = nan([numel(reg_lab) numel(cond_lab) numel(SBJs)]);
for s = 1:numel(SBJs)
    % Load RL Model
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
    
    % Average within condition
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        model(:,cond_ix,s) = nanmean(tmp.model(full_cond_idx==cond_ix,:),1);
    end
end

%% Compute Group Averages
% Variance is very small, so plot standard deviation instead of standard error of the mean
plot_means = nanmean(model,3);
plot_stds  = nan([numel(reg_lab) numel(cond_lab)]);
% plot_sems  = nan([numel(reg_lab) numel(cond_lab)]);
for reg_ix = 1:numel(reg_lab)
    plot_stds(reg_ix,:) = nanstd(model(reg_ix,:,:),[],3);
%     plot_sems(reg_ix,:) = nanstd(model(reg_ix,:,:),[],3)./sqrt(numel(SBJs))';
end

%% Plot Predictors
fig_name = [SBJ_id '_' model_id '_predictions'];
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','OuterPosition',[0 0 0.5 0.5]);
ax = gca; hold on;

% Plot Easy and Hard as separate lines
cond_x = 1:numel(cond_lab);
ez_ix = cond_x(ez_idx);
hd_ix = cond_x(~ez_idx);

reg_lines = gobjects(size(reg_lab));
for reg_ix = 1:numel(reg_lab)
    % Plot group mean
    reg_lines(reg_ix) = plot(ez_ix+plt.x_fudge*(reg_ix-1),plot_means(reg_ix,ez_ix),'Color',reg_colors{reg_ix},...
        'LineStyle',reg_styles{reg_ix},'LineWidth',plt.width);%,'Marker',reg_mrkrs{reg_ix});
    plot(hd_ix+plt.x_fudge*(reg_ix-1),plot_means(reg_ix,hd_ix),'Color',reg_colors{reg_ix},...
        'LineStyle',reg_styles{reg_ix},'LineWidth',plt.width);%,'Marker',reg_mrkrs{reg_ix});
    
    % Plot variance
    errorbar(ez_ix+plt.x_fudge*(reg_ix-1),plot_means(reg_ix,ez_ix),plot_stds(reg_ix,ez_ix),...%plot_sems
        'Color',reg_colors{reg_ix},'LineStyle',reg_styles{reg_ix},'LineWidth',plt.width);
    errorbar(hd_ix+plt.x_fudge*(reg_ix-1),plot_means(reg_ix,hd_ix),plot_stds(reg_ix,hd_ix),...%plot_sems
        'Color',reg_colors{reg_ix},'LineStyle',reg_styles{reg_ix},'LineWidth',plt.width);
end

% Plot parameters
set(gca,'XTick',1:numel(cond_lab));
set(gca,'XTickLabels',cond_names);
% xtickangle(plt.tick_angle);
xlim([0 numel(cond_lab)+1]);
if any(strcmp(st.model_lab,{'VML','SML'}))
    ylim([-1.1 1.1]);
    plt.leg_loc = 'southwest';
end
% yticks([0 1]);
% set(gca,'YTickLabels',{'Low','High'})

if strcmp(st.model_lab,'RSVPE')
    model_str = 'Reward Valence vs. Value vs. Prediction Error';
elseif any(strcmp(st.model_lab,{'VML','SML'}))
    model_str = 'Outcome Features';
elseif any(strcmp(st.model_lab,{'ERPEsL','RPEsL'}))
    model_str = 'Reward Prediction Error (RPE) Features';
else
    model_str = model_id;
end
title(model_str,'Interpreter','none');
legend(reg_lines,reg_names,'Location',plt.leg_loc);
set(gca,'FontSize',16);

%% Save Figure
if save_fig
    % Create figure directory
    fig_dir = [root_dir 'PRJ_Error_eeg/results/model_predictions/' plt_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
