function SBJ04e_ERP_plot_RL_elec_comparison_R2_ts(SBJ_id,an_ids,stat_id,plt_id,save_fig,varargin)
%% Plots Adjusted R2 model fits across different electrodes (same model) for time series
%   Option: Select original or adjusted R2
%   Only for single channel right now...
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   an_ids [cell array] - string IDs of the analysis parameters to compare
%   stat_id [str] - ID of the stats parameters to use
%   plt_id [str] - ID of the plotting parameters to use
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
%       r2_version [str] - {'Adjusted' or 'Ordinary'} version of R2
%           default: 'Adjusted'
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
        elseif strcmp(varargin{v},'r2_version')
            % 'Ordinary' or 'Adjusted' (default)
            r2_version = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('r2_version','var'); r2_version = 'Adjusted'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Analysis and Plotting Parameters
if numel(an_ids)>2; error('only ready for 2 elecs now'); end
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get Plotting Parameters
load([root_dir 'PRJ_Error_EEG/data/' SBJs{1} '/04_proc/' SBJs{1} '_' an_ids{1} '.mat'],'roi');
cfgs = []; cfgs.latency = st.stat_lim;
st_roi = ft_selectdata(cfgs, roi);
st_time_vec = st_roi.time{1};
ch_list = st_roi.label;
an_styles = {'-',':'};

%% Load Models
lmes = cell([numel(an_ids) numel(st_time_vec)]);
for an_ix = 1:numel(an_ids)
    tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_ids{an_ix} '.mat']);
    lmes(an_ix,:) = tmp.lme;
end

%% Plot Model Comparisons
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' strjoin(an_ids,'-') '/' stat_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Get R2 and mean
    r2s = NaN([numel(an_ids) numel(st_time_vec)]);
    mean_r2 = NaN(size(an_ids));
    an_leg  = cell(size(an_ids));
    for an_ix = 1:numel(an_ids)
        for t_ix = 1:numel(st_time_vec)
            r2s(an_ix,t_ix) = lmes{an_ix,t_ix}.Rsquared.(r2_version);
        end
        mean_r2(an_ix) = nanmean(r2s(an_ix,:));
        an_leg{an_ix} = [an_ids{an_ix} ' (mean=' num2str(mean_r2(an_ix),'%.3f') ')'];
    end
    
    %% Create plot
    if strcmp(r2_version,'Ordinary')
        fig_name = [SBJ_id '_RL_R2ord_comparison_' stat_id];
    else
        fig_name = [SBJ_id '_RL_R2adj_comparison_' stat_id];
    end
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
    
    %% Plot R2
    ax = gca; hold on;
    
    % Plot R2 per electrode
    main_lines = gobjects(size(an_ids));
    for an_ix = 1:numel(an_ids)
        main_lines(an_ix) = line(st_time_vec, r2s(an_ix,:),...
            'Color','k','LineWidth',2,'LineStyle',an_styles{an_ix});
    end
    ylims = ylim;
    
    % Axes and Labels
    ax.YLabel.String = [r2_version ' R2'];
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    title([ch_list{ch_ix} ' (n=' num2str(numel(SBJs)) ')']);
    if plt.legend
        legend(main_lines,an_leg,'Location',plt.legend_loc,'Interpreter','none');
    end
    set(gca,'FontSize',16);
    
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
