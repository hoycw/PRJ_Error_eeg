function fn_plot_ERP_stack(SBJ, proc_id, plt_id, data, fig_vis, save_fig, path)
%% Plot ERP with stacked single trials and mean for ICA components
%   Currently only plots stimulus-locked, full trial data
% INPUTS:
%   SBJ [str] - name of the SBJ
%   proc_id [str] - name of the preprocessing pipeline parameters (e.g., 'egg_full_ft')
%   plt_id [str] - name of plotting parameter set
%   data [FT struct] - data to plot (cut into trials)
%   fig_vis [0/1] - if a data_browser view of the time course of the ICA
%   save_fig [0/1] - binary flag to save figure
%   path [str] - path name to directory in which "path/plot/SBJ_chlab_ERP_stack.png" will save

%% Data Preparation
% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
plt_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_cmd);

%% Trim data to plotting time
cfg = [];
cfg.latency = plt.plt_lim;
data = ft_selectdata(cfg,data);

%% Compute ERP
cfg = [];
cfg.channel = 'all';
erps = ft_timelockanalysis(cfg, data);

%% Event processing
% Load behavioral params
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

if ~strcmp(proc.event_type,'S')
    % !!! improve logic here to look for compatibility between plt events and what's available based on trial_lim_s in proc_vars
    error('mismatch in events in plt and proc_vars!');
end

evnt_ix = zeros([3 1]);
% Stimulus
evnt_ix(1) = find(data.time{1}==0);
% Response
evnt_ix(2) = find(data.time{1}==prdm_vars.target);
% Feedback
evnt_ix(3) = find(data.time{1}==prdm_vars.target+prdm_vars.fb_delay);

%% Plot Data
% keep manual screen position - better in dual monitor settings

for ch_ix = 1:numel(data.label)
    fig_name = [SBJ '_' data.label{ch_ix}];
    figure('units','normalized','Name', fig_name, 'Visible', fig_vis);
    subplot('Position', [0.1, 0.35, 0.8, 0.6]);
    
    % Get data matrix
    plot_data = NaN([numel(data.trial) numel(data.time{1})]);
    for t_ix = 1:numel(data.trial)
        plot_data(t_ix, :) = squeeze(data.trial{t_ix}(ch_ix,:));
    end
    
    % Get color limits
    %clims = NaN([1 2]);
    %clims(1) = prctile(plot_data(:),5);
    %clims(2) = prctile(plot_data(:),95);
    
    % Plot single trial stack
    imagesc(plot_data);
    set(gca,'YDir','normal');
    
    % Plot events
    for e_ix = 1:numel(plt.evnt_type)
        line([evnt_ix(e_ix) evnt_ix(e_ix)],ylim,...
            'LineWidth',plt.evnt_width(e_ix),'Color',plt.evnt_color{e_ix},'LineStyle',plt.evnt_style{e_ix});
    end
    
    % Plot Labels and Parameters
    ax = gca;
    ax.YLabel.String = 'Trials';
    ax.XLim          = [0,size(data.trial{1},2)];
    ax.XTick         = 0:plt.x_step_sz*data.fsample:size(data.trial{1},2);
    ax.XTickLabel    = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    title(data.label{ch_ix});
    if plt.legend
        legend(plt.evnt_type,'Location',plt.legend_loc);
    end
    colorbar;
    %caxis(clims);
    
    % Plot ERP
    subplot('Position', [0.1, 0.1, 0.70, 0.15]);
    plot(erps.time,erps.avg(ch_ix,:),'k');
    xlim(plt.plt_lim);
    
    % Save figure
    if save_fig
        comp_stack_fname = [path 'plot/' SBJ data.label{ch_ix} '_ERP_stack.png'];
        saveas(gcf,comp_stack_fname);
    end
end

end
