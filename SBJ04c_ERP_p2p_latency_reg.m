function SBJ04c_ERP_p2p_latency_reg(SBJ_id,proc_id,an_id,pk_stat_id,SBJ_norm,save_fig,varargin)
% Run and plot Linear Mixed-Effects Model on FRN peak latency using condition-averaged RL model predictors
%   Must run SBJ04c_ERP_grp_stats_LME_P2P first to obtain peak data
%   Only for single channel
% COMPUTATIONS:
%   Load single-trial design matrix (model regressors) and average within condtion
%   Load peak times identified in peak-to-peak FRN LME analysis
%   LME multiple regression predicting peak latency using model regressors
%   Scatter plot of latencies with simple linear fit for visualization
%       Plotting Option: normalize latencies within SBJ by subtracting mean latency across conditions
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   pk_stat_id [str] - ID of the peak-to-peak stats analysis to provide peak data
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure (and prints outcome statistics)

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
if ~exist('fig_vis','var');        fig_vis = 'on'; end
if ~exist('fig_ftype','var');      fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Analysis and Plotting Parameters
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' pk_stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.measure,'p2p') || ~strcmp(st.an_style,'lme'); error('run only for p2p LME analyses!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, reg_names, ~, ~]  = fn_regressor_label_styles(st.model_lab);
[cond_lab, cond_names, cond_colors, ~, cond_markers] = fn_condition_label_styles(st.trial_cond{1});

%% Compute mean regressor per condition
cond_reg_mean = nan([numel(cond_lab) numel(SBJs) numel(reg_lab)]);
plot_reg_mean = nan([numel(cond_lab) numel(SBJs) numel(reg_lab)]);
for s = 1:numel(SBJs)
    % Load data
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
    
    % Z-score SBJ model regressors
    plot_model = tmp.model;
    sbj_model = NaN(size(tmp.model));
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            sbj_model(:,reg_ix) = ...
                (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
        end
    else
        sbj_model = tmp.model;
    end
    
    % Compute mean regressor per condition
    cond_idx = fn_condition_index(cond_lab, bhv);
    for reg_ix = 1:numel(reg_lab)
        for cond_ix = 1:numel(cond_lab)
            cond_reg_mean(cond_ix,s,reg_ix) = nanmean(sbj_model(cond_idx==cond_ix,reg_ix));
            % Keep non-normalized regressors for plotting
            plot_reg_mean(cond_ix,s,reg_ix) = nanmean(plot_model(cond_idx==cond_ix,reg_ix));
        end
    end
    
    clear tmp
end

%% Load Peak Times
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' pk_stat_id '_' an_id '.mat']);
if numel(ch_list)>1; error('only plotting for 1 channel in this script!'); end
if st.pk_sign(2)~=-1; error('second peak is not negative!'); end

% Take second (negative) peak; assume 1 channel
pk_data = squeeze(pk_times(:,:,1,2));

% Normalize peak latencies within SBJ
if SBJ_norm
    pk_data = pk_data-mean(pk_data,1);
    norm_str = '_SBJnorm';
else
    norm_str = '';
end

%% Comptue latency-predictor regression
% Build Model Table
tbl = table;
tbl.latency = pk_data(:);
if any(isnan(pk_data(:))); warning([num2str(sum(isnan(pk_data(:)))) ' NaNs in LME peak data!']);end
for reg_ix = 1:numel(reg_lab)
    tbl.(reg_lab{reg_ix}) = reshape(cond_reg_mean(:,:,reg_ix),numel(cond_lab)*numel(SBJs),1);
end

% Create Model Formula
reg_formula = strjoin(reg_lab,' + ');
formula = ['latency ~ ' reg_formula];

glm = fitglm(tbl,formula);
pvals = glm.Coefficients.pValue(2:end); % Skip intercept

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(pvals);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

%% Plot latencies and simple regression
% Figure parameters
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' pk_stat_id '/pk_lat_reg/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
fig_name = [SBJ_id '_pk_lat_reg_' pk_stat_id '_' an_id '_' ch_list{1} norm_str];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.8],'Visible',fig_vis);

for reg_ix = 1:numel(reg_lab)
    subplot(2,2,reg_ix);
    ax = gca; hold on;
    
    % Plot Latency Data
    scats = gobjects(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        scats(cond_ix) = scatter(pk_data(cond_ix,:),plot_reg_mean(cond_ix,:,reg_ix),50,...
            cond_colors{cond_ix},cond_markers{cond_ix});
    end
    
    % Compute simple linear fit for visualization
    tmp_reg = reshape(plot_reg_mean(:,:,reg_ix),numel(cond_lab)*numel(SBJs),1);
    tmp_pk  = pk_data(:);
    lin_fit = polyfit(tmp_pk(~isnan(tmp_pk)),tmp_reg(~isnan(tmp_pk)),1);
    reg_x_fudge = 0.001;
    reg_x_step  = 0.001;
    reg_x = min(pk_data(:))-reg_x_fudge:reg_x_step:max(pk_data(:))+reg_x_fudge;
    reg_y = lin_fit(1)*reg_x + lin_fit(2);
    
    % Plot linear fit
    reg_line = plot(reg_x,reg_y,'k');
    
    % Scale Y axis to range of model regressors
    switch reg_lab{reg_ix}
        case 'EV'
            ax.YLim = [-1 1];
        case 'sRPE'
            ax.YLim = [-2 2];
        case 'uRPE'
            ax.YLim = [0 2];
        case 'Lik'
            ax.YLim = [0 1];
        otherwise
            error(['Unknown reg_lab: ' reg_lab{reg_ix}]);
    end
    
    % Plotting parameters
    ax.YLabel.String = reg_names{reg_ix};
    if SBJ_norm
        line([0 0],ylim,'Color','k','LineStyle','--');
        ax.XLabel.String = 'Normalized FRN Peak Latency (s)';
        ax.XLim = [-0.06 0.06];
        ax.XTick = [-0.05:0.025:0.05];
    else
        ax.XLabel.String = 'FRN Peak Latency (s)';
        ax.XLim = st.pk_lim(2,:);
    end
    ax.Title.String = [ch_list{1} ' ' reg_lab{reg_ix} ': beta=' num2str(glm.Coefficients.Estimate(reg_ix+1),'%.03f')...
        ' (p=' num2str(qvals(reg_ix),'%.03f') ')'];
    legend(scats,cond_names,'Location','best');
    set(ax,'FontSize',16');
    
    % Print multiple regression stats
    fprintf('%s p = %f\n',reg_lab{reg_ix},qvals(reg_ix));
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
