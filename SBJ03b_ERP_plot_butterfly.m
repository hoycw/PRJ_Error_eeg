function SBJ03b_ERP_plot_butterfly(SBJ,conditions,proc_id,an_id,plt_id,save_fig,varargin)
%% Plot single trial data (butterfly) for single SBJ, one condition per subplot
% INPUTS:
%   SBJ [str] - ID of subject to run
%   conditions [str] - group of condition labels to segregate trials
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
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
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

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
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load data
load([SBJ_vars.dirs.proc SBJ '_' an_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

% Select conditions (and trials)
[cond_lab, cond_names, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(cond_lab, bhv);

% Get trials for plotting
trials = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(roi.label) numel(cond_trial_ix) numel(roi.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
    end
end

%% Get event timing for plotting
[evnt_times] = fn_get_evnt_times(an.event_type,plt.evnt_lab,'prdm_vars',prdm_vars);

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' conditions '/but_' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(roi.label)
    %% Compute plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(roi.time{1})]);
    sems  = NaN([numel(cond_lab) numel(roi.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(ch_ix,:,:),2));
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(ch_ix,:,:),[],2))./sqrt(size(trials{cond_ix},2))';
    end
    
    %% Create plot
    fig_name = [SBJ '_' conditions '_' an_id '_' roi.label{ch_ix} '_but'];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 1],'Visible',fig_vis);
    
    %% Butterfly Plots
    axes = gobjects([numel(cond_lab)+1 1]);
    for cond_ix = 1:numel(cond_lab)
        axes(cond_ix) = subplot(numel(cond_lab)+1, 1, cond_ix); hold on;
        
        % Plot individual trials per condition
        plot(roi.time{1},squeeze(trials{cond_ix}(ch_ix,:,:)),...
            'Color',cond_colors{cond_ix},'LineWidth',plt.butterfly_width,...
            'LineStyle',cond_styles{cond_ix});
        
        % Plot Events
        for evnt_ix = 1:numel(plt.evnt_lab)
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        % Axes and Labels
        axes(cond_ix).YLabel.String = 'uV';
        axes(cond_ix).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
        axes(cond_ix).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
        axes(cond_ix).Title.String  = [roi.label{ch_ix} ': ' cond_lab{cond_ix} ' (n=' ...
            num2str(size(trials{cond_ix},2)) ')'];
    end
    
    %% Plot Means (and variance)
    axes(end) = subplot(numel(cond_lab)+1, 1, numel(cond_lab)+1); hold on;
    
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(plt.evnt_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(roi.time{1}, means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_lines(numel(cond_lab)+evnt_ix) = line(...
            [evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end
    leg_lab = [cond_lab plt.evnt_lab];
    
    % Axes and Labels
    axes(end).YLabel.String = 'uV';
    axes(end).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(end).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(end).XLabel.String = 'Time (s)';
    axes(end).Title.String  = roi.label{ch_ix};
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    
    %% Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
