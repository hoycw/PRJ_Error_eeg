function SBJ04e_ERP_plot_RL_model_comparison_point(SBJ_id,an_id,stat_ids,null_id,plt_id,save_fig,varargin)
% Plots AIC model fits across different RL models for point estimates
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
    if strcmp(st.measure,'ts'); error('this script is for point estimates!');end
    
    if st_ix>1
        fnames = fieldnames(sts{st_ix});
        for f_ix = 1:numel(fnames)
            if ~any(strcmp(fnames{f_ix},{'model_lab','z_reg'}))
                if ischar(sts{st_ix}.(fnames{f_ix}))
                    if ~strcmp(sts{1}.(fnames{f_ix}), sts{st_ix}.(fnames{f_ix}))
                        error(['st.' fnames{f_ix} ' not the same!']);
                    end
                elseif isnumeric(sts{st_ix}.(fnames{f_ix}))
                    if any(sts{1}.(fnames{f_ix}) ~= sts{st_ix}.(fnames{f_ix}))
                        error(['st.' fnames{f_ix} ' not the same!']);
                    end
                elseif iscell(sts{st_ix}.(fnames{f_ix}))
                    if numel(sts{st_ix}.(fnames{f_ix}))>1; error('not ready for cells n > 1'); end
                    if any(sts{1}.(fnames{f_ix}){1} ~= sts{st_ix}.(fnames{f_ix}){1})
                        error(['st.' fnames{f_ix} ' not the same!']);
                    end
                else
                    error(['Unknown class(' fnames{f_ix} ')']);
                end
            end
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
ch_list = roi.label;
if numel(ch_list)>1; error('only for 1 channel now...'); end
st_colors = distinguishable_colors(numel(stat_ids)+plot_null);

%% Load Models
% Load real models
aics = NaN(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_ids{st_ix} '_' an_id '.mat']);
    aics(st_ix) = tmp.lme{1}.ModelCriterion.AIC;
    if st_ix==1
        if strcmp(sts{st_ix}.measure,'mean')
            st_lim = sts{st_ix}.stat_lim + tmp.reg_pk_time;
            measure_str = [sts{1}.measure '(' num2str(st_lim(1)) '-' num2str(st_lim(2)) ')'];
        elseif strcmp(sts{st_ix}.measure,'p2p')
            error('write p2p measure_str');
        else
            error(['Unknown st.measure: ' sts{st_ix}.measure]);
        end
    end
end

% Load null model
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' null_id '_' an_id '.mat']);
null_aic = tmp.lme{1}.ModelCriterion.AIC;

%% Plot Model Comparisons
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' strjoin(stat_ids,'-') '/' null_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Compute plotting data
% Compute relative likelihoods
aic_min = min(aics);
rel_lik = nan(size(stat_ids));
rl_leg  = cell(size(stat_ids));
model_labs = cell(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    rel_lik(st_ix) = exp((aic_min-aics(st_ix))/2);
    model_labs{st_ix} = sts{st_ix}.model_lab;
    rl_leg{st_ix} = ['RL=' num2str(rel_lik(st_ix),'%.2f')];
end

% Compute for null model
rel_lik_null = exp((aic_min-null_aic)/2);
null_rl_leg = ['RL=' num2str(rel_lik_null,'%.2f')];

%% Create plot
fig_name = [SBJ_id '_RL_AIC_comparison_' an_id];
if plot_null
    fig_name = [fig_name '_null'];
end
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots

%% Plot AIC
ax = gca; hold on;

% Create plotting variables
if plot_null
    plot_ids = [model_labs {'SBJonly'}];
    plot_aics = [aics null_aic];
    plot_RLs = [rl_leg null_rl_leg];
    leg = [rl_leg null_rl_leg];
else
    plot_ids = model_labs;
    plot_aics = aics;
    plot_RLs = rl_leg;
    leg = rl_leg;
end

% Get YLim
rl_fudge = plt.sig_yfudge/2;
y_range = max(plot_aics)-min(plot_aics);

% Plot AIC
b = bar(1:numel(plot_ids),plot_aics,'BarWidth',plt.bar_width,'FaceColor','flat');
% for st_ix = 1:numel(plot_ids)
%     b.CData(st_ix,:) = st_colors(st_ix,:);
% end

% Add Relative Likelihoods
for st_ix = 1:numel(plot_ids)
    text(st_ix,plot_aics(st_ix)+y_range*rl_fudge,plot_RLs{st_ix},'FontSize',16,'HorizontalAlignment','center');
end

% Axes and Labels
ax.YLim          = [min(plot_aics)-y_range*plt.sig_yfudge max(plot_aics)+y_range*plt.sig_yfudge];
ax.YLabel.String = 'AIC';
ax.XLim          = [0 numel(plot_ids)+1];
ax.XTick         = 1:numel(plot_ids);
ax.XTickLabel    = plot_ids;
ax.XLabel.String = 'Models';
title([ch_list{1} ' ' measure_str]);
% if plt.legend
%     legend(b,leg,'Location','best','Interpreter','none');%plt.legend_loc
% end
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
