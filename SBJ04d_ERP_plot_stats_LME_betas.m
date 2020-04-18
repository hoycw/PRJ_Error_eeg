function SBJ04d_ERP_plot_stats_LME_betas(SBJ_id,an_id,stat_id,plt_id,save_fig,varargin)
% Plots beta weights per regressor for mean window LME analyses (e.g., FRN)
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
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Analysis and Plotting Parameters
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.measure,'mean'); error('run only for mean window LME analyses!'); end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = load_SBJ_file(SBJ_id);

% Select Conditions of Interest
[reg_lab, reg_names, reg_colors, ~]  = fn_regressor_label_styles(st.model_lab);

%% Load Stats
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat']);
if numel(ch_list)>1; error('only plotting for 1 channel in this script!'); end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' stat_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Compute plotting data
% Obtain Model Parameters
if strcmp(st.model_lab,'SBJonly')
    plot_betas = lme{1}.Coefficients.Estimate;
else
    plot_betas = lme{1}.Coefficients.Estimate(2:end);
end
r2 = lme{1}.Rsquared.Adjusted;

% Find significant regressors
sig_reg    = false(size(reg_lab));
for reg_ix = 1:numel(reg_lab)
    if any(qvals(reg_ix,:) <= st.alpha)
        sig_reg(reg_ix) = true;
    end
end

%% Create plot
fig_name = [SBJ_id '_' stat_id '_' ch_list{1}];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
ax = gca; hold on;

% Find significance marker distance
sig_y = max(abs(plot_betas))*plt.sig_yfudge;

% Plot Betas
b = bar(1:numel(reg_lab),plot_betas,'BarWidth',plt.bar_width,'FaceColor','flat');
for reg_ix = 1:numel(reg_lab)
    b.CData(reg_ix,:) = reg_colors{reg_ix};
    
    % Plot Significance
    if sig_reg(reg_ix)
        if qvals(reg_ix)<=plt.sig_cut3
            x_pos = [-plt.sig_xfudge 0 plt.sig_xfudge];
        elseif qvals(reg_ix)<=plt.sig_cut2
            x_pos = [-plt.sig_xfudge plt.sig_xfudge]/2;
        else
            x_pos = 0;
        end
        y_pos = plot_betas(reg_ix)+sig_y*sign(plot_betas(reg_ix));
        scatter(x_pos+reg_ix,repmat(y_pos,size(x_pos)),...
            plt.sig_sz,plt.sig_color,plt.sig_mrkr);
    end
end

% Axes and Labels
%     ax.YLim          = ylims; %!!! change for plt.sigType=line
ax.YLabel.String = 'Beta Weight';
ax.XLim          = [0 numel(reg_lab)+1];
ax.XTick         = 1:numel(reg_lab);
ax.XTickLabel    = reg_names;
ax.XLabel.String = 'Regressors';
ax.Title.String  = [ch_list{1} ': R2 = ' num2str(r2,'%.3f') ' (n = ' num2str(numel(SBJs)) ')'];
if plt.legend
    legend(b,reg_lab,'Location',plt.legend_loc);
end
ylims = ylim;
set(gca,'FontSize',16);
ax.YLim = ylims;

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