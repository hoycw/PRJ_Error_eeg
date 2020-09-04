function SBJ03b_ERP_plot_stack(SBJ,conditions,proc_id,an_id,plt_id,save_fig,varargin)
%% Plot ERPs and single trial stack for single SBJ
%   Intended to be a QA tool (e.g., plot win vs. loss without easy/hard)
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

% Load data, behavior, model
load([SBJ_vars.dirs.proc,SBJ,'_',an_id,'.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

% Select conditions (and trials)
[cond_lab, ~, cond_colors, cond_styles, cond_markers] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(cond_lab, bhv);

% Split trials by win/loss for plotting
stack_cond = zeros([sum(cond_idx~=0) 1]);
stack = NaN([numel(roi.label) sum(cond_idx~=0) numel(roi.time{1})]);
stack_cnt_ix = 0;
for cond_ix = 1:numel(cond_lab)
    % Extract trials per condition
    stack_trial_ix = find(cond_idx==cond_ix);
    for trl_ix = 1:numel(stack_trial_ix)
        stack_cnt_ix = stack_cnt_ix + 1;
        stack(:,stack_cnt_ix,:) = roi.trial{stack_trial_ix(trl_ix)};
        stack_cond(stack_cnt_ix) = cond_ix;
    end
end

%% Get event timing for plotting
[evnt_times] = fn_get_evnt_times(an.event_type,plt.evnt_lab,'prdm_vars',prdm_vars);

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' conditions '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(roi.label)
    %% Compute plotting data    
    % Compile stack
    means = NaN([numel(cond_lab) numel(roi.time{1})]);
    sems  = NaN([numel(cond_lab) numel(roi.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(stack(ch_ix,stack_cond==cond_ix,:),2));
        sems(cond_ix,:) = squeeze(std(stack(ch_ix,stack_cond==cond_ix,:),[],2))./sqrt(sum(stack_cond==cond_ix))';
    end
    
    %% Create plot
    fig_name = [SBJ '_' conditions '_' an_id '_' roi.label{ch_ix}];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 1],'Visible',fig_vis);
    
    %% Plot Trial Stack
    axes = gobjects([2 1]);
    axes(1) = subplot(3,1,[1 2]); hold on;
    
    % Plot single trial stack
    imagesc(roi.time{1},1:size(stack,2),squeeze(stack(ch_ix,:,:)));
    set(gca,'YDir','Normal');
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        if any(strcmp(plt.evnt_lab{evnt_ix},{'F','Fon'}))
            for trl_ix = 1:size(stack,2)
                scatter(evnt_times(evnt_ix), trl_ix,...
                    'MarkerFaceColor',[cond_colors{stack_cond(trl_ix)}],'MarkerEdgeColor','none',...
                    'Marker','o');%cond_markers{stack_cond(trl_ix)});
            end
        else
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],[1 size(stack,2)],...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_style);
        end
    end
    
    % Axes and labels
    axes(1).YLabel.String = 'Trials';
    axes(1).YLim = [1 size(stack,2)];
    axes(1).XLim = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(1).XLabel.String = 'Time (s)';
    title(roi.label{ch_ix});
%     if plt.legend
%         legend(an.event_type,'Location',plt.legend_loc);
%     end
    colorbar('Location','northoutside');
    
    %% Plot ERPs
    axes(2) = subplot(3,1,3); hold on;
    
    % Plot individual trials per condition under ERPs
    if plt.butterfly
        for cond_ix = 1:numel(cond_lab)
            plot(roi.time{1},squeeze(trials{cond_ix}(ch_ix,:,:)),...
                'Color',cond_colors{cond_ix},'LineWidth',plt.butterfly_width,...
                'LineStyle',cond_styles{cond_ix});
        end
    end
    
    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(roi.time{1}, means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_lines(numel(cond_lab)+evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end
    leg_lab = [cond_lab an.event_type];
        
    % Axes and Labels
    axes(2).YLabel.String = 'uV';
    axes(2).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(2).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(2).XLabel.String = 'Time (s)';
    axes(2).Title.String  = [roi.label{ch_ix} ': ' conditions];
    if plt.legend
        legend(main_lines,leg_lab{:},'Location','best');%plt.legend_loc);
    end
    
    %% Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
