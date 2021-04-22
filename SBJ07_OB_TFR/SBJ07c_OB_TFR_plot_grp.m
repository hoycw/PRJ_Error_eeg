function SBJ07c_OB_TFR_plot_grp(SBJ_id,conditions,proc_id,an_id,plt_id,save_fig,varargin)
%% Plot TFRs of power data per condition averaged across the group
% INPUTS:
%   SBJ_id [str] - ID of group of subjects to plot
%   conditions [str] - group of condition labels to segregate trials
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
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

%% Load Results
if ~strcmp(proc_id,'odd_full_ft'); error('SBJ07 only for oddball task!'); end
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select conditions (and trials)
[cond_lab, cond_names, ~, ~, ~] = fn_condition_label_styles(conditions);

% Load data
tfr_all = cell([numel(cond_lab) numel(SBJs)]);
for s = 1:numel(SBJs)
    fprintf('-------------------- Processing %s ------------------------\n',SBJs{s});
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    tmp = load([SBJ_vars.dirs.proc SBJs{s} '_' proc_id '_' an_id '.mat']);
    
    % Average across trials within condition/subject
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

%% Get event timing for plotting
if ~strcmp(an.event_type,'S') || ~strcmp(plt.evnt_lab,'S'); error('only S for oddball!'); end
evnt_times = 0;

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' an_id '/' conditions '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr_avg{1}.label)
    %% Create plot
    fig_name = [SBJ_id '_' conditions '_' an_id '_' tfr_avg{1}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.7],'Visible',fig_vis);
    
    % Get color lims per condition
    clim = zeros([numel(cond_lab) 2]);
    for cond_ix = 1:numel(cond_lab)
        clim(cond_ix,:) = [min(tfr_avg{cond_ix}.powspctrm(:)) max(tfr_avg{cond_ix}.powspctrm(:))];
    end
    tick_ix = 1:3:numel(tfr_avg{1}.freq);
    yticklab = cell(size(tick_ix));
    for f = 1:numel(tick_ix)
        yticklab{f} = num2str(tfr_avg{1}.freq(tick_ix(f)),'%.1f');
    end
    
    % Condition Plots
    %   Originally used ft_singleplotTFR, but that wasn't as flexible
    %   Switched to manual plotting with imagesc
    [n_rowcol, ~] = fn_num_subplots(numel(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        subplot(n_rowcol(1),n_rowcol(2),cond_ix);
        
        % Plot TFR
        imagesc(tfr_avg{cond_ix}.time, 1:numel(tfr_avg{cond_ix}.freq),...
            squeeze(tfr_avg{cond_ix}.powspctrm(ch_ix,:,:)));%,[min(clim(:,1)) max(clim(:,2))]);
        set(gca,'YDir','normal');
%         set(gca,'YTick',1:3:numel(tfr_avg{cond_ix}.freq));
%         set(gca,'YTickLabels',yticklab);
        caxis([-max(abs(clim(:))) max(abs(clim(:)))]);
        %ft_singleplotTFR(cfgplt, tfr_avg{cond_ix});
        
        % Plot Events
        for evnt_ix = 1:numel(plt.evnt_lab)
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        % Axes and parameters
        title([tfr_avg{cond_ix}.label{ch_ix} ': ' cond_names{cond_ix} '(n=' num2str(numel(SBJs)) ')']);
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
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
