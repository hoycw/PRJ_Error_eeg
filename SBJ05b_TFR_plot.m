function SBJ05b_TFR_plot(SBJ,conditions,proc_id ,an_id, save_fig,varargin)
%% Plot ERPs for single SBJ
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
elseif strcmp(varargin{v},'write_report')
write_report = varargin{v+1};
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
%plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
%eval(plt_vars_cmd);

% Load data
load([SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',an_id,'.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

% Select conditions (and trials)
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(cond_lab, bhv);


%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr.label)
    %% Create plot
        fig_name = [SBJ '_' conditions '_' an_id '_' ...
            tfr.label{ch_ix}];
        figure('Name',fig_name,'units','normalized',...
       'outerposition',[0 0 1 0.5],'Visible',fig_vis);
        % Get trials for plotting
        tfr_conds = cell(size(cond_lab));
        for cond_ix = 1:numel(cond_lab)
            cond_trial_ix = find(cond_idx==cond_ix);
            tfr_conds{cond_ix} = tfr;
            tfr_conds{cond_ix}.powspctrm = tfr_conds{cond_ix}.powspctrm(cond_trial_ix, ch_ix, :,:);
            tfr_conds{cond_ix}.trialinfo = tfr_conds{cond_ix}.trialinfo(cond_trial_ix);
            tfr_conds{cond_ix}.cumtapcnt = tfr_conds{cond_ix}.cumtapcnt(cond_trial_ix);
        end
        cfg_s = [];
        cfg_s.channel = tfr.label(ch_ix);
        % Condition Plots
        for cond_ix = 1:length(cond_lab)
            subplot(1,numel(cond_lab)+1,cond_ix);   %!!! assumes only 2 conditions
            cfgraw = [];
            cfgraw.baseline = an.bsln_lim; %should be in sec
            cfgraw.trials = 'all';
            cfgraw.xlim = an.trial_lim_s;
            ft_singleplotTFR(cfgraw, tfr_conds{cond_ix});
            ax = gca;
            title(cond_lab{cond_ix});
        end
end
%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',an_id, conditions, '_conds.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'tfr_conds');
% Save figure
if save_fig
fig_filename = [fig_dir fig_name '.png'];
fprintf('Saving %s\n',fig_filename);
saveas(gcf,fig_filename);
%eval(['export_fig ' fig_filename]);
end

