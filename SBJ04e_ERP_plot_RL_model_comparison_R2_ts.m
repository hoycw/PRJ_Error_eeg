function SBJ04e_ERP_plot_RL_model_comparison_R2_ts(SBJ_id,an_id,stat_ids,null_id,plt_id,save_fig,varargin)
% Plots R2 model fits across different RL models for time series
%   Option: Select original or adjusted R2
%   Option: If null_id is not empty (''), subtract off R2 for that stat_id
%   Only for single channel right now...
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   an_id [str] - ID of the analysis parameters to use
%   stat_ids [cell array] - string IDs of the stats parameters to compare
%   null_id [str] - ID of the SBJonly baseline model to compare
%   plt_id [str] - ID of the plotting parameters to use
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
%       rm_null [0/1] - binary flag to subtract null model R2 as baseline
%           default: 1
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
        elseif strcmp(varargin{v},'rm_null')
            rm_null = varargin{v+1};
        elseif strcmp(varargin{v},'r2_version')
            % 'Ordinary' or 'Adjusted' (default)
            r2_version = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');    fig_vis = 'on'; end
if ~exist('fig_ftype','var');  fig_ftype = 'png'; end
if ~exist('rm_null','var');    rm_null = 1; end
if ~exist('r2_version','var'); r2_version = 'Adjusted'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Analysis and Plotting Parameters
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Load stat parameters and check compatibility
sts = cell(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_ids{st_ix} '_vars.m'];
    eval(stat_vars_cmd);
    sts{st_ix} = st;
    
    % Check alignment of time windows and measurements
    if st_ix>1
        if any(sts{1}.stat_lim ~= sts{st_ix}.stat_lim)
            error('st.stat_lim not aligned!');
        end
        if ~strcmp(sts{1}.measure, sts{st_ix}.measure)
            error('st.measure not the same!');
        end
    end
    clear st stat_vars_cmd
end

% Load SBJonly null model
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' null_id '_vars.m'];
eval(stat_vars_cmd);
null_st = st;

% Check compatibility of null model
if any(sts{1}.stat_lim ~= null_st.stat_lim)
    error('st.stat_lim not aligned!');
end
if ~strcmp(sts{1}.measure, null_st.measure)
    error('st.measure not the same!');
end
clear st stat_vars_cmd

% Get Plotting Parameters
load([root_dir 'PRJ_Error_EEG/data/' SBJs{1} '/04_proc/' SBJs{1} '_' an_id '.mat'],'roi');
cfgs = []; cfgs.latency = sts{1}.stat_lim;
st_roi = ft_selectdata(cfgs, roi);
st_time_vec = st_roi.time{1};
ch_list = st_roi.label;
st_colors = distinguishable_colors(numel(stat_ids));

%% Load Models
% Load real models
lmes = cell([numel(stat_ids) numel(st_time_vec)]);
for st_ix = 1:numel(stat_ids)
    tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_ids{st_ix} '_' an_id '.mat']);
    lmes(st_ix,:) = tmp.lme;
end

% Load null model
null_r2 = zeros(size(st_time_vec));
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' null_id '_' an_id '.mat']);
for t_ix = 1:numel(st_time_vec)
    null_r2(t_ix) = tmp.lme{t_ix}.Rsquared.(r2_version);
end

%% Plot Model Comparisons
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/model_comparisons/' strjoin(stat_ids,'-') '/' null_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Compute R2 and mean
    r2s     = NaN([numel(stat_ids) numel(st_time_vec)]);
    mean_r2 = NaN(size(stat_ids));
    st_leg  = cell(size(stat_ids));
    for st_ix = 1:numel(stat_ids)
        for t_ix = 1:numel(st_time_vec)
            if rm_null
                r2s(st_ix,t_ix) = lmes{st_ix,t_ix}.Rsquared.(r2_version)-null_r2(t_ix);
            else
                r2s(st_ix,t_ix) = lmes{st_ix,t_ix}.Rsquared.(r2_version);
            end
        end
        mean_r2(st_ix) = nanmean(r2s(st_ix,:));
        st_leg{st_ix} = [stat_ids{st_ix} ' (mean=' num2str(mean_r2(st_ix),'%.3f') ')'];
    end
    
    % Compute for null model
    mean_r2_null = nanmean(null_r2);
    null_leg = [null_id ' (mean=' num2str(mean_r2_null,'%.3f') ')'];
    
    %% Create plot
    if strcmp(r2_version,'Ordinary')
        fig_name = [SBJ_id '_RL_R2ord_comparison_' an_id];
    else
        fig_name = [SBJ_id '_RL_R2adj_comparison_' an_id];
    end
    if rm_null
        fig_name = [fig_name '_rmnull'];
    end
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
    
    %% Plot R2
    ax = gca; hold on;
    
    % Plot R2 per model
    main_lines = gobjects(size(stat_ids));
    for st_ix = 1:numel(stat_ids)
        main_lines(st_ix) = line(st_time_vec, r2s(st_ix,:),...
            'Color',st_colors(st_ix,:),'LineWidth',2);
    end
    
    % Plot null model R2
    if ~rm_null
        main_lines(end+1) = plot(st_time_vec, null_r2,...
            'Color', [0.4 0.4 0.4], 'LineStyle', '--');
        st_leg = [st_leg null_leg];
    end
    ylims = ylim;
    
    % Axes and Labels
    ax.YLabel.String = [r2_version ' R2'];
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    if rm_null
        title([ch_list{ch_ix} ' (n=' num2str(numel(SBJs)) '); SBJ null removed']);
    else
        title([ch_list{ch_ix} ' (n=' num2str(numel(SBJs)) ')']);
    end
    if plt.legend
        legend(main_lines,st_leg,'Location',plt.legend_loc,'Interpreter','none');
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
