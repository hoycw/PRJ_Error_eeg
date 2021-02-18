function SBJ03c_ERP_plot_diff_grp(SBJ_id,conditions,proc_id,an_id,plt_id,save_fig,varargin)
%% Plot ERPs and difference waves averaged across group
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
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
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select conditions (and trials)
[cond_lab, cond_names, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);
[reg_lab, ~, reg_colors, ~]  = fn_regressor_label_styles('ERPEsL');

% Create contrast:
[diff_lab, diff_names, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles(conditions);

% Load example data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m'];
eval(SBJ_vars_cmd);
tmp = load([root_dir 'PRJ_Error_eeg/data/',SBJs{1},'/04_proc/',SBJs{1},'_',an_id,'.mat'],'roi');
time_vec = tmp.roi.time{1};
ch_list  = tmp.roi.label;

% Load all data
means      = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
diff_waves = nan([numel(diff_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    load([SBJ_vars.dirs.proc,SBJs{s},'_',an_id,'.mat']);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    
    % Select Conditions of Interest
    cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for t_ix = 1:numel(cond_trial_ix)
            trials(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
        end
        means(cond_ix,s,:,:) = mean(trials,2);
    end
    
    % Compute SBJ difference waves
    for pair_ix = 1:numel(diff_lab)
        diff_waves(pair_ix,s,:,:) = means(diff_pairs(pair_ix,1),s,:,:)-means(diff_pairs(pair_ix,2),s,:,:);
    end
    
    clear tmp SBJ_vars bhv roi
end

%% Get event timing for plotting
if strcmp(an.event_type,'S')
    error('add loading of prdm_vars to plot relative to stim!');
end
[evnt_times] = fn_get_evnt_times(an.event_type,plt.evnt_lab);

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' conditions '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Compute ERP means and variance
    erp_means = NaN([numel(cond_lab) numel(time_vec)]);
    erp_sems  = NaN([numel(cond_lab) numel(time_vec)]);
    for cond_ix = 1:numel(cond_lab)
        erp_means(cond_ix,:) = squeeze(mean(means(cond_ix,:,ch_ix,:),2));
        erp_sems(cond_ix,:) = squeeze(std(means(cond_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    
    % Compute difference wave means and variance
    diff_means = NaN([numel(diff_lab) numel(time_vec)]);
    diff_sems  = NaN([numel(diff_lab) numel(time_vec)]);
    for pair_ix = 1:numel(diff_lab)
        diff_means(pair_ix,:) = squeeze(mean(diff_waves(pair_ix,:,ch_ix,:),2));
        diff_sems(pair_ix,:) = squeeze(std(diff_waves(pair_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    
    %% Print summary stats
    fprintf('-------- Difference Wave Properties ----------\n');
    for pair_ix = 1:numel(diff_lab)
        % Report max
        [~, max_ix] = max(diff_means(pair_ix,:));
        fprintf('\t(%d) %s [max, time] = [%.2f at %.3f]\n',pair_ix,diff_names{pair_ix},...
            diff_means(pair_ix,max_ix),time_vec(max_ix));
        
        % Report min
        [~, min_ix] = min(diff_means(pair_ix,:));
        fprintf('\t(%d) %s [min, time] = [%.2f at %.3f]\n',pair_ix,diff_names{pair_ix},...
            diff_means(pair_ix,min_ix),time_vec(min_ix));
    end
    
    % Report max/min for mean across condition pairs
    mean_diff = mean(diff_means,1);
    [~, max_ix] = max(mean_diff);
    [~, min_ix] = min(mean_diff);
    fprintf('\t%s Overall (max, time) = [%.2f at %.3f]\n',conditions,mean_diff(max_ix),time_vec(max_ix));
    fprintf('\t%s Overall (min, time) = [%.2f at %.3f]\n',conditions,mean_diff(min_ix),time_vec(min_ix));

    %% Create ERP plot
    fig_name = [SBJ_id '_' conditions '_' an_id '_' ch_list{ch_ix}];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.7 1.0],'Visible',fig_vis);
    subplot(2,1,1); ax = gca; hold on;
    
    % Plot ERP Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(time_vec, erp_means(cond_ix,:), erp_sems(cond_ix,:),...
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
    leg_lab = [cond_names plt.evnt_lab];
    
    % Axes and Labels
    ax.YLabel.String = 'uV';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    ax.Title.String  = [ch_list{ch_ix} ' ERPs (n=' num2str(numel(SBJs)) ')'];
    set(ax,'FontSize',16');
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    
    %% Plot Difference Waves
    subplot(2,1,2); ax = gca; hold on;
    
    % Set color to match RL model and ylim for goodall Fz/Pz comparison
    switch conditions
        case 'RewP'
            ylims = [-4 10];
        case 'Pos-Neg'
            ylims = [-5 20];
            mean_color = reg_colors{strcmp(reg_lab,'sRPE')};
        case 'Large-Small'
            ylims = [-10 10];
            mean_color = reg_colors{strcmp(reg_lab,'uRPE')};
        case 'Unlik-Lik'
            ylims = [-6 10];
            mean_color = reg_colors{strcmp(reg_lab,'Lik')};
    end
    if ~strcmp(SBJ_id,'goodall'); clear ylims; end
    
    % Plot Difference Wave Grand Average (and variance)
    ebars = cell(size(diff_lab));
    main_lines = gobjects([numel(diff_lab)+numel(an.event_type) 1]);
    for pair_ix = 1:numel(diff_lab)
        ebars{pair_ix} = shadedErrorBar(time_vec, diff_means(pair_ix,:), diff_sems(pair_ix,:),...
            'lineProps',{'Color',diff_colors{pair_ix},'LineWidth',plt.mean_width,...
            'LineStyle',diff_styles{pair_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(pair_ix) = ebars{pair_ix}.mainLine;
    end
    
    % Plot Mean Difference Wave
    if numel(diff_lab)>1
        mn_diff_plot = 1;
        main_lines(numel(diff_lab)+mn_diff_plot) = line(time_vec, mean_diff, ...
            'Color', mean_color, 'LineWidth', 3);
        leg_lab = [diff_names 'Mean'];
    else
        mn_diff_plot = 0;
        leg_lab = diff_names;
    end
    if ~exist('ylims','var'); ylims = ylim; end
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_lines(numel(diff_lab)+mn_diff_plot+evnt_ix) = line(...
            [evnt_times(evnt_ix) evnt_times(evnt_ix)],ylims,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end
    leg_lab = [leg_lab plt.evnt_lab];
    
    % Add zero line
    plot(xlim,[0 0],'Color',[0 0 0],'LineWidth',0.2);
    
    % Axes and Labels
    ax.YLabel.String = 'uV';
    ax.YLim          = ylims;
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    ax.Title.String  = [ch_list{ch_ix} ' ' conditions ' (n=' num2str(numel(SBJs)) ')'];
    set(ax,'FontSize',16');
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
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

end
