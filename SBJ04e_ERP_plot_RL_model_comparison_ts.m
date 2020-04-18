function SBJ04e_ERP_plot_RL_model_comparison_ts(SBJ_id,an_id,stat_ids,null_id,plt_id,save_fig,varargin)
% Plots AIC model fits across different RL models for time series
%   Only for single channel right now...
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
        elseif strcmp(varargin{v},'plot_null')
            plot_null = varargin{v+1};
%         elseif strcmp(varargin{v},'r2_version')
%             % 'Ordinary' or 'Adjusted' (default)
%             r2_version = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('plot_null','var'); plot_null = 1; end
% if ~exist('r2_version','var'); r2_version = 'Adjusted'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Analysis and Plotting Parameters
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = load_SBJ_file(SBJ_id);

sts = cell(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_ids{st_ix} '_vars.m'];
    eval(stat_vars_cmd);
    sts{st_ix} = st;
    if ~strcmp(st.measure,'ts'); error('this script is for time series!');end
    
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

if any(sts{1}.stat_lim ~= null_st.stat_lim)
    error('null_st.stat_lim not aligned!');
end
if ~strcmp(sts{1}.measure, null_st.measure)
    error('null_st.measure not the same!');
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
null_aic = zeros(size(st_time_vec));
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' null_id '_' an_id '.mat']);
for t_ix = 1:numel(st_time_vec)
    null_aic(t_ix) = tmp.lme{t_ix}.ModelCriterion.AIC;
end

%% Plot Model Comparisons
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' strjoin(stat_ids,'-') '/' null_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Collect AIC
    aics = NaN([numel(stat_ids) numel(st_time_vec)]);
    mean_aic = NaN(size(stat_ids));
    for st_ix = 1:numel(stat_ids)
        for t_ix = 1:numel(st_time_vec)
            aics(st_ix,t_ix) = lmes{st_ix,t_ix}.ModelCriterion.AIC;
        end
        mean_aic(st_ix) = nanmean(aics(st_ix,:));
    end
    
    % Compute relative likelihoods
    aic_min = min(mean_aic);
    rel_lik = nan(size(stat_ids));
    st_leg  = cell(size(stat_ids));
    for st_ix = 1:numel(stat_ids)
        rel_lik(st_ix) = exp((aic_min-mean_aic(st_ix))/2);
        st_leg{st_ix} = [stat_ids{st_ix} ' (mean=' num2str(round(mean_aic(st_ix)))...
            '; RL=' num2str(rel_lik(st_ix),'%.2f') ')'];
    end
    
    % Compute for null model
    mean_aic_null = nanmean(null_aic);
    rel_lik_null = exp((aic_min-mean_aic_null)/2);
    null_leg = [null_id ' (mean=' num2str(round(mean_aic_null)) ...
            '; RL=' num2str(rel_lik_null,'%.2f') ')'];
    
    %% Create plot
    fig_name = [SBJ_id '_RL_AIC_comparison_' an_id];
    if plot_null
        fig_name = [fig_name '_null'];
    end
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
    
    %% Plot AIC
    ax = gca; hold on;
    
    % Plot Means (and variance)
    main_lines = gobjects(size(stat_ids));
    for st_ix = 1:numel(stat_ids)
        main_lines(st_ix) = line(st_time_vec, aics(st_ix,:),...
            'Color',st_colors(st_ix,:),'LineWidth',2);
    end
    if plot_null
        main_lines(end+1) = plot(st_time_vec, null_aic,...
            'Color', [0.4 0.4 0.4], 'LineStyle', '--');
        st_leg = [st_leg null_leg];
    end
    ylims = ylim;
    
    % Axes and Labels
    ax.YLabel.String = 'AIC';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    title([ch_list{ch_ix} ' (n=' num2str(numel(SBJs)) ')']);
    if plt.legend
        legend(main_lines,st_leg,'Location','best','Interpreter','none');%plt.legend_loc
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
