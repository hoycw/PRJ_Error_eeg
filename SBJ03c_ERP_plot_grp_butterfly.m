function SBJ03c_ERP_plot_grp_butterfly(SBJs,conditions,proc_id,an_id,plt_id,save_fig,varargin)
%% Plot ERPs for group, with butterfly plots for SBJ means per condition
% INPUTS:
%   conditions [str] - group of condition labels to segregate trials

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
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
%         elseif strcmp(varargin{v},'write_report')
%             write_report = varargin{v+1};
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

[cond_lab, cond_names, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);

% Load example data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m'];
eval(SBJ_vars_cmd);
tmp = load([root_dir 'PRJ_Error_eeg/data/',SBJs{1},'/04_proc/',SBJs{1},'_',an_id,'.mat'],'roi');
time_vec = tmp.roi.time{1};
ch_list  = tmp.roi.label;

% Load all data
means = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
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
    
    clear tmp SBJ_vars bhv roi
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' conditions '/' an_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Compute means and variance
    plt_means = NaN([numel(cond_lab) numel(time_vec)]);
    sems  = NaN([numel(cond_lab) numel(time_vec)]);
    for cond_ix = 1:numel(cond_lab)
        plt_means(cond_ix,:) = squeeze(mean(means(cond_ix,:,ch_ix,:),2));
        sems(cond_ix,:) = squeeze(std(means(cond_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    
    %% Create plot
    fig_name = ['GRP_' conditions '_' an_id '_' ch_list{ch_ix} '_but'];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 1],'Visible',fig_vis);   %this size is for single plots
%     [plot_rc,~] = fn_num_subplots(numel(w2.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    axes = gobjects([numel(cond_lab)+1 1]);
    subplot(numel(cond_lab)+1,1,1);
    axes(1) = gca; hold on;
    
    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(time_vec, plt_means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    % Plot Extra Features (events, significance)
    for evnt_ix = 1:numel(an.event_type)
        main_lines(numel(cond_lab)+evnt_ix) = line([0 0],ylim,...
            'LineWidth',plt.evnt_width(evnt_ix),'Color',plt.evnt_color{evnt_ix},...
            'LineStyle',plt.evnt_style{evnt_ix});
    end
    leg_lab = [cond_lab an.event_type];
    
    % Axes and Labels
    axes(1).YLabel.String = 'uV';
    axes(1).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(1).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(1).XLabel.String = 'Time (s)';
    axes(1).Title.String  = [ch_list{ch_ix} ' (n=' num2str(numel(SBJs)) ')'];
    set(axes(1),'FontSize',16');
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    
    %% Plot Butterfly for Each Condition
    for cond_ix = 1:numel(cond_lab)
        subplot(numel(cond_lab)+1,1,cond_ix+1);
        axes(cond_ix+1) = gca; hold on;
        
        % Plot SBJ Means
        plot(time_vec, squeeze(means(cond_ix,:,ch_ix,:)),'Color',cond_colors{cond_ix},...
            'LineWidth',plt.butterfly_width);
        
        % Plots GRP Mean
        plot(time_vec, plt_means(cond_ix,:), 'Color', 'k', 'LineWidth', plt.mean_width);
        
        % Plot Extra Features (events, significance)
        for evnt_ix = 1:numel(an.event_type)
            line([0 0],ylim,'LineWidth',plt.evnt_width(evnt_ix),...
                'Color',plt.evnt_color{evnt_ix},'LineStyle',plt.evnt_style{evnt_ix});
        end
        % leg_lab = [cond_lab an.event_type];
        
        % Axes and Labels
        axes(cond_ix+1).YLabel.String = 'uV';
        axes(cond_ix+1).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
        axes(cond_ix+1).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
        axes(cond_ix+1).XLabel.String = 'Time (s)';
        axes(cond_ix+1).Title.String  = [ch_list{ch_ix} ': ' cond_lab{cond_ix}];
        set(axes(cond_ix+1),'FontSize',16');
%         if plt.legend
%             legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
%         end
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
