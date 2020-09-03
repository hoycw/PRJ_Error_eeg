function SBJ04d_ERP_plot_stats_LME_RL_topo_reg(SBJ_id,an_id,stat_id,plt_id,save_fig,varargin)
% Plots group RL beta topographies per regressor for single ERP analysis (time point)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
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

%% Analysis and Plotting Parameters
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select Conditions of Interest
[reg_lab, ~, reg_colors, reg_styles]  = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, cond_colors, cond_styles, ~] = fn_condition_label_styles(st.trial_cond{1});

% Check for window compatibility
if strcmp(st.measure,'ts')
    error('Single topo plot must have single metric within window!');
end

%% Load Stats
% Check stat analysis was run on correct SBJs
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'SBJs');
if ~all(strcmp(SBJs,tmp.SBJs))
    fprintf(2,'Loaded SBJs: %s\n',strjoin(tmp.SBJs,', '));
    error('Not all SBJs match input SBJ list!');
end

% Load results
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'lme','qvals','ch_list','reg_pk_time');
if numel(ch_list)<64; error('Cannot plot opo wihtout full cap!'); end

% Get beta values and color limits
clim  = zeros([1 2]);
betas = nan(size(qvals));
r2    = zeros(size(ch_list));
for ch_ix = 1:numel(ch_list)
    for reg_ix = 1:numel(reg_lab)
        betas(reg_ix,ch_ix) = lme{ch_ix}.Coefficients.Estimate(reg_ix+1);
        clim = [min([clim(1) betas(reg_ix,ch_ix)]) max([clim(2) betas(reg_ix,ch_ix)])];
    end
    r2(ch_ix) = lme{ch_ix}.Rsquared.Adjusted;
end

%% Create dummy dataset for plotting
% This will initialize electrode labels, neighbor structure, etc.
load([root_dir 'PRJ_Error_eeg/data/',SBJs{1},'/04_proc/',SBJs{1},'_',an_id,'.mat'],'roi');
topo = {};
topo.label  = roi.label;
topo.time   = reg_pk_time;
topo.dimord = 'chan_time';

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' stat_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create plot
fig_name = [SBJ_id '_' stat_id '_' an_id];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);

% Initialize plotting parameters
cfgp = [];
cfgp.xlim            = [reg_pk_time reg_pk_time];
cfgp.zlim            = clim;
cfgp.layout          = 'biosemi64.lay';
cfgp.colorbar        = 'yes';
cfgp.comment         = 'no';
cfgp.highlight       = 'on';
cfgp.highlightsymbol = '*';
cfgp.maskparameter   = 'mask';

% Create a subplot for each regressor
axes = gobjects([numel(reg_lab)+1 1]);
[num_rc,~] = fn_num_subplots(numel(reg_lab)+1);
for reg_ix = 1:numel(reg_lab)
    subplot(num_rc(1),num_rc(2),reg_ix);
    axes(reg_ix) = gca; hold on;
    
    % Insert data and mark significant electrodes
    topo.avg  = betas(reg_ix,:)';
    ch_ix = 1:numel(topo.label);
    cfgp.highlightchannel = ch_ix(qvals(reg_ix,:)'<=st.alpha);
    topo.mask = ones(size(topo.avg))*0.1;
    topo.mask(qvals(reg_ix,:)'<=st.alpha) = 1;
    % cfgp.zlim = [min(betas(reg_ix,:)) max(betas(reg_ix,:))];
    
    % Plot Beta Topos
    ft_topoplotER(cfgp, topo);
    title(reg_lab{reg_ix});
    axis tight        
end

% Plot R2 Topo
subplot(num_rc(1),num_rc(2),numel(reg_lab)+1);
axes(end) = gca; hold on;

topo.avg  = r2;
topo.mask = ones(size(topo.avg));
cfgp.zlim = [min(r2) max(r2)];
ft_topoplotER(cfgp, topo);
title(['Adjusted R2 (' num2str(numel(SBJs)) ')']);
axis tight

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
