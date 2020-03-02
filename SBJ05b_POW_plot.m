function SBJ05b_POW_plot(SBJ,conditions,proc_id,an_id,plt_id,save_fig,varargin)
%% Plot band-limited power for single SBJ
% INPUTS:
%   conditions [str] - group of condition labels to segregate trials

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
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
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);
load([SBJ_vars.dirs.proc SBJ '_' proc_id '_' an_id '.mat']);
if numel(tfr.freq)>1; error('TFR is not averaged over frequencies!'); end

% Select conditions (and trials)
[cond_lab, ~, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(cond_lab, bhv);

%% Get event timing for plotting
evnt_times = zeros(size(plt.evnt_lab));
if strcmp(an.event_type,'S')
    for evnt_ix = 1:numel(plt.evnt_lab)
        switch plt.evnt_lab{evnt_ix}
            case 'S'
                evnt_times(evnt_ix) = 0;
            case 'R'
                evnt_times(evnt_ix) = prdm_vars.target;
            case {'F','Fon'}
                evnt_times(evnt_ix) = prdm_vars.target+prdm_vars.fb_delay;
            case 'Foff'
                evnt_times(evnt_ix) = prdm_vars.target+prdm_vars.fb_delay+prdm_vars.fb;
            otherwise
                error(['Unknown event type in plt: ' plt.evnt_lab{evnt_ix}]);
        end
    end
elseif strcmp(an.event_type,'F')
    evnt_times(1) = 0;
else
    error('Unknown an.event_type');
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' an_id '/' conditions '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr.label)
    %% Compute plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(tfr.time)]);
    sems  = NaN([numel(cond_lab) numel(tfr.time)]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(tfr.powspctrm(cond_idx==cond_ix,ch_ix,:,:),1));
        sems(cond_ix,:)  = squeeze(std(tfr.powspctrm(cond_idx==cond_ix,ch_ix,:,:),[],1))./sqrt(sum(cond_idx==cond_ix))';
    end
    
    %% Create plot
    fig_name = [SBJ '_' conditions '_' an_id '_' tfr.label{ch_ix}];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
%     [plot_rc,~] = fn_num_subplots(numel(tfr.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    ax = gca; hold on;
    
    % Plot individual trials per condition
    if plt.butterfly
        error('butterfly not done yet!');
%         for cond_ix = 1:numel(cond_lab)
%             plot(tfr.time,squeeze(trials{cond_ix}(ch_ix,:,:)),...
%                 'Color',cond_colors{cond_ix},'LineWidth',plt.butterfly_width,...
%                 'LineStyle',cond_styles{cond_ix});
%         end
    end
    
    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(plt.evnt_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(tfr.time, means(cond_ix,:), sems(cond_ix,:),...
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
    ax.YLabel.String = 'uV';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    ax.Title.String  = tfr.label{ch_ix};
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
