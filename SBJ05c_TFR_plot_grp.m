function SBJ05c_TFR_plot_grp(SBJs,conditions,proc_id,an_id,save_fig,varargin)
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
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Select conditions (and trials)
[grp_lab, ~, ~] = fn_group_label_styles(conditions);
[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);
% if ~strcmp(st.model_lab,{'DifOut','Out'}); error('not ready for surprise trials!'); end
grp_cond_lab = cell(size(grp_lab));
for grp_ix = 1:numel(grp_lab)
    [grp_cond_lab{grp_ix}, ~, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
end

% Load data
tfr_all = cell([numel(cond_lab) numel(SBJs)]);
for s = 1:numel(SBJs)
    fprintf('-------------------- Processing %s ------------------------\n',SBJs{s});
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    tmp = load([SBJ_vars.dirs.proc,SBJs{s},'_',an_id,'.mat']);
    for cond_ix = 1:numel(cond_lab)
        cfgs = [];
        cfgs.trials = find(fn_condition_index(cond_lab(cond_ix), bhv));
        cfgs.avgovertrials = 'yes';
        tfr_all{cond_ix, s} = ft_freqdescriptives(cfgs, tmp.tfr);
    end
    clear bhv tfr SBJ_vars
end

% Average across SBJs
tfr_avg = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    tfr_avg{cond_ix} = ft_freqgrandaverage([], tfr_all{cond_ix,:});
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr_avg{1}.label)
    %% Create plot
    fig_name = ['GRP_' conditions '_' an_id '_' tfr_avg{1}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Get color lims per condition
    clim = zeros([numel(cond_lab) 2]);
    for cond_ix = 1:numel(cond_lab)
        clim(cond_ix,:) = [min(tfr_all{cond_ix}.powspctrm(:)) max(tfr_all{cond_ix}.powspctrm(:))];
    end
    tick_ix = 1:3:numel(tfr_all{1}.freq);
    yticklab = cell(size(tick_ix));
    for f = 1:numel(tick_ix)
        yticklab{f} = num2str(tfr_all{1}.freq(tick_ix(f)),'%.1f');
    end
    
    % Condition Plots
    %cfgplt = []; cfgplt.zlim = clim;
    for cond_ix = 1:length(cond_lab)
        subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix);
        imagesc(tfr_all{cond_ix}.time, 1:numel(tfr_all{cond_ix}.freq), squeeze(tfr_all{cond_ix}.powspctrm(ch_ix,:,:)),[min(clim(:,1)) max(clim(:,2))]);
        set(gca,'YDir','normal');
        set(gca,'YTick',1:3:numel(tfr_all{cond_ix}.freq));
        set(gca,'YTickLabels',yticklab);
        %ft_singleplotTFR(cfgplt, tfr_avg{cond_ix});
        title([tfr_avg{cond_ix}.label{ch_ix} ': ' cond_lab{cond_ix}]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar;
        set(gca,'FontSize',16);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.png'];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end

end
