function SBJ05c_TFR_plot_grp(SBJs,conditions, an_id, save_fig,varargin)
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
for sbj_ix = 1:numel(SBJs)
    SBJ_vars_cmd{sbj_ix} = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{sbj_ix} '_vars.m'];
    eval(SBJ_vars_cmd{sbj_ix});
    SBJ_vars_all{sbj_ix} = SBJ_vars;
end

an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);

% Load data
tfr_conds_all = cell(size(SBJs));
for sbj_ix = 1:numel(SBJs)
    tmp = load([SBJ_vars_all{sbj_ix}.dirs.SBJ,'04_proc/',SBJs{sbj_ix},'_',an_id,conditions, '_conds.mat']); tfr_conds_all{sbj_ix} = tmp.tfr_conds;
end
cfg_avg = [];
for cond_ix = 1:numel(cond_lab)
    for sbj_ix = 1: numel(SBJs)
        averages{sbj_ix}= tfr_conds_all{sbj_ix}{cond_ix};
        tmp_tfr{sbj_ix} = ft_freqdescriptives(cfg_avg, averages{sbj_ix});
    end
    tfr_avg{cond_ix} = ft_freqgrandaverage(cfg_avg, tmp_tfr{1:numel(SBJs)});
end
%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/GRP/' an_id '/'];
if ~exist(fig_dir,'dir')
mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr_avg{1}.label)
    %% Create plot
        fig_name = ['GRP_' conditions '_' an_id '_' ...
            tfr_avg{1}.label{ch_ix}];
        figure('Name',fig_name,'units','normalized',...
       'outerposition',[0 0 1 0.5],'Visible',fig_vis);
        % Get trials for plotting
        cfg_s = [];
        cfg_s.channel = tfr_avg{1}.label{ch_ix};
        % Condition Plots
        for cond_ix = 1:numel(cond_lab)
            subplot(1,numel(cond_lab)+1,cond_ix);   %!!! assumes only 2 conditions
            cfgraw = [];
            cfgraw.baseline = an.bsln_lim; %should be in sec
            cfgraw.xlim = an.trial_lim_s;
            ft_singleplotTFR(cfgraw, tfr_avg{cond_ix});
            ax = gca;
            title(cond_lab{cond_ix});
        end
end
% Save figure
if save_fig
fig_filename = [fig_dir fig_name '.png'];
fprintf('Saving %s\n',fig_filename);
saveas(gcf,fig_filename);
%eval(['export_fig ' fig_filename]);
end

