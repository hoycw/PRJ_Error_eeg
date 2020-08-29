function SBJ04b_BHV_RL_model_plot_grp(SBJ_id,proc_id,stat_id, varargin)
%% Plot all SBJ RL model fits overlapping
%   Rereferenced to the midpoint tolerance
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       save_fig [0/1] - binary flag to save figure; default = 1
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig') && ischar(varargin{v+1})
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');     fig_vis = 'on'; end
if ~exist('fig_ftype','var');   fig_ftype = 'png'; end
if ~exist('save_fig','var');    save_fig = 1; end

%% Load Parameters
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Load SBJ list
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
model_id = [st.model_lab '_' st.trial_cond{1}];
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

% Initialize Plotting Variables
sig_step = 0.001;               % tolerance step size for plotting model fit
sig_x    = [0:sig_step:0.4];
sig_y    = nan([numel(sig_x) numel(SBJs)]);

%% Load Data
for s = 1: numel(SBJs)
    % Load Subject Specific Data
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.proc SBJ '_model_' model_id '.mat']);
    load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'],'bhv');

    % Select Conditions of Interest
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    bhv_fields = fieldnames(bhv);
    orig_n_trials = numel(bhv.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
            bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
        end
    end
    n_trials = numel(bhv.trl_n);
    fprintf('\t%s: Loaded %d trials, kept %d for modeling...\n',SBJ,orig_n_trials,n_trials);
    
    % Recalculate Model Fit
    sig_y(:,s) = betas(1) + (sig_x * betas(2));
    sig_y(:,s) = 1 ./ (1+exp(-sig_y(:,s)));
end

%% Plot All SBJ Level Sigmoids
fig_name = ['GRP_BHV_acc_' model_id '_pWin'];
fig = figure('Name',fig_name,'Visible',fig_vis);
hold on;

% Plot model fits
for s = 1:numel(SBJs)
    line(sig_x,sig_y(:,s),'Color', 'k');
end

% Figure Parameters
title('GRP Level Accuracy vs. Tolerance');
xlabel('Tolerance (s)');
ylabel('Accuracy');
set(gca,'YLim',[0 1]);
set(gca,'FontSize',14);

%% Save Figures
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/model_fits/' model_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    % Ensure vector graphics if saving
    if any(strcmp(fig_ftype,{'svg','eps'}))
        set(gcf, 'Renderer', 'painters');
    end
    saveas(fig,fig_fname);
end

end
