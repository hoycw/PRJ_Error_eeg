function fn_plot_ERP_stack(SBJ, proc_id, plt_id, data, fig_vis, save_fig)
%% Plot ERP with stacked single trials for ICA components

%% Data Preparation
% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Trim data to plotting time
cfg = [];
cfg.latency = plt_vars.plt_lim;
data = ft_selectdata(cfg,data);

%% Compute ERP
cfg = [];
cfg.channel = 'all';
erps = ft_timelockanalysis(cfg, data);

%% Event processing
% Load behavioral params
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

if ~strcmp(proc_vars.event_type,'S')
    % !!! improve logic here to look for compatibility between plt_vars events and what's availabel based on trial_lim_s in proc_vars
    error('mismatch in events in plt_vars and proc_vars!');
end

evnt_ix = zeros([3 1]);
% Stimulus
evnt_ix(1) = find(data.time{1}==0);
% Response
evnt_ix(2) = find(data.time{1}==prdm_vars.target);
% Feedback
evnt_ix(3) = find(data.time{1}==prdm_vars.target+prdm_vars.fb_delay);

%% Import behavioral data
% %   Total_Trial,Block,Condition,Hit,RT,Timestamp,Tolerance,Trial,Score,ITI,ITI type
% 
% events = ft_read_event(SBJ_vars.dirs.raw_filename);
% eventvalues = extractfield(events, 'value');
% eventsamples = extractfield(events, 'sample');
% if numel(eventvalues) ~= numel(eventsamples)
%    difference = numel(eventsamples)-numel(eventvalues);
%    eventsamples = eventsamples(difference+1:end);
% end
% stimulus_ix = find(eventvalues == 1);
% feedback_ix = find(eventvalues == 2);
% stimulus_onset = eventsamples(stimulus_ix);
% feedback_onset = eventsamples(feedback_ix);
% stim_onset = find(data.time{1,1} == 0);
% stim_onset_array(1:numel(stimulus_onset)) = stim_onset;
% [~,fb_onset] = min(abs(1.8-data.time{1,1}));
% fb_onset_array(1:numel(feedback_onset)) = fb_onset;
% [~,one_sec] = min(abs(1 - data.time{1,1}));
% one_sec_array(1:numel(feedback_onset)) = one_sec;
% %for x = 1: numel(feedback_onset);
%     %trial_sample = feedback_onset(x)/numel(data.time{1,x}(1,:))

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
    % SHEILA: we actually weren't using this below, so you computed it and
    % left it... however, maybe it's nice ot see outliers for data
    % cleaning, so try with and without and decide which is better.
    %clims = NaN([1 2]);
    %clims(1) = prctile(plot_data(:),5);
    %clims(2) = prctile(plot_data(:),95);

    
    % Plot single trial stack
    imagesc(plot_data);
    set(gca,'YDir','normal');
    
    % Plot events
    for e_ix = 1:numel(plt_vars.evnt_type)
        line([evnt_ix(e_ix) evnt_ix(e_ix)],ylim,...
            'LineWidth',plt_vars.evnt_width(e_ix),'Color',plt_vars.evnt_color{e_ix},'LineStyle',plt_vars.evnt_style{e_ix});
    end
    
    ax = gca;
%     ax.XTick = [1, 257, 513, 769, 1025, 1281, 1537, 1793, 2049, 2305, 2561, 2817, 3073];
%     ax.XTickLabel = [-.25, 0, 0.25, .5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75];
    ax.YLabel.String = 'Trials';
    ax.XLim          = [0,size(data.trial{1},2)];
    ax.XTick         = 0:plt_vars.x_step_sz*data.fsample:size(data.trial{1},2);
    ax.XTickLabel    = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    title(data.label{ch_ix});
    if plt_vars.legend
        legend(plt_vars.evnt_type,'Location',plt_vars.legend_loc);
    end
    colorbar;
    %caxis(clims);
    
    % Plot ERP
    subplot('Position', [0.1, 0.1, 0.70, 0.15]);
    plot(erps.time,erps.avg(ch_ix,:),'k');
    xlim(plt_vars.plt_lim);
    
    % Save figure
    if save_fig
        comp_stack_fname = [SBJ_vars.dirs.proc_stack SBJ data.label{ch_ix} '_ERP_stack.png'];  % SHEILA: add a sub directory for ERPstack 
        saveas(gcf,comp_stack_fname);
    end
end

end