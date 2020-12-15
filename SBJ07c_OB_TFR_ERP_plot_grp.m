function SBJ07c_OB_TFR_ERP_plot_grp(SBJ_id,conditions,proc_id,tfr_an_id,erp_an_id,plt_id,save_fig,varargin)
%% Plot TFRs of power data with ERP overlap per condition averaged across the group
% INPUTS:
%   SBJ_id [str] - ID of group of subjects to plot
%   conditions [str] - group of condition labels to segregate trials
%   proc_id [str] - ID of preprocessing pipeline
%   tfr_an_id [str] - ID of the TFR analysis parameters to use
%   erp_an_id [str] - ID of the ERP analysis parameters to use
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
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' erp_an_id '_vars.m'];
eval(an_vars_cmd);
erp_an = an;
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' tfr_an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~strcmp(an.event_type,erp_an.event_type); error('TFR and erp event mismatch!'); end
if ~all(an.trial_lim_s==erp_an.trial_lim_s); error('TFR and erp trial_lim_s mismatch!'); end
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
    load([SBJ_vars.dirs.proc SBJs{s} '_' erp_an_id '.mat'],'roi');
    load([SBJ_vars.dirs.proc SBJs{s} '_' proc_id '_' tfr_an_id '.mat'],'tfr');
    
    % Initialize matrices
    if s==1
        roi_time_vec = roi.time{1};
        means    = nan([numel(cond_lab) numel(SBJs) numel(tfr.label) numel(roi_time_vec)]);
    end
    
    % Average across trials within condition/subject
    cfgs = []; cfgs.avgovertrials = 'yes';
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(fn_condition_index(cond_lab(cond_ix), bhv));
        cfgs.trials = cond_trial_ix;
        tfr_all{cond_ix, s} = ft_freqdescriptives(cfgs, tfr);
        
        % ERP Computation
        trials = nan([numel(tfr.label) numel(cond_trial_ix) numel(roi_time_vec)]);
        for t_ix = 1:numel(cond_trial_ix)
            trials(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
        end
        means(cond_ix,s,:,:) = mean(trials,2);
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
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' tfr_an_id '/' conditions '/' erp_an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr_avg{1}.label)
    %% Compute ERP plotting data    
    % Compute means and variance
    plt_means = NaN([numel(cond_lab) numel(roi_time_vec)]);
    sems  = NaN([numel(cond_lab) numel(roi_time_vec)]);
    for cond_ix = 1:numel(cond_lab)
        plt_means(cond_ix,:) = squeeze(mean(means(cond_ix,:,ch_ix,:),2));
        sems(cond_ix,:) = squeeze(std(means(cond_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    erp_ylim = [min(min(plt_means-sems)) max(max(plt_means+sems))];
    
    %% Create plot
    fig_name = [SBJ_id '_' conditions '_' tfr_an_id '_' erp_an_id '_' tfr_avg{1}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);%[0 0 0.8 0.8]
    
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
    %cfgplt = []; cfgplt.zlim = clim;
    [n_rowcol, ~] = fn_num_subplots(numel(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        subplot(n_rowcol(1),n_rowcol(2),cond_ix);
        
        % Plot TFR
        yyaxis left
        imagesc(tfr_avg{cond_ix}.time, tfr_avg{cond_ix}.freq, squeeze(tfr_avg{cond_ix}.powspctrm(ch_ix,:,:)));%1:numel(tfr_all{cond_ix}.freq)
        set(gca,'YDir','normal');
%         set(gca,'YTick',1:3:numel(tfr_all{cond_ix}.freq));
%         set(gca,'YTickLabels',yticklab);
        caxis([-max(abs(clim(:))) max(abs(clim(:)))]);
        colorbar('northoutside');
        ylabel('Frequency (Hz)');

        % Plot Events
        for evnt_ix = 1:numel(plt.evnt_lab)
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        % Plot Means (and variance)
        yyaxis right
        shadedErrorBar(roi_time_vec, plt_means(cond_ix,:), sems(cond_ix,:),...
                'lineProps',{'Color','k','LineWidth',2,...
                'LineStyle','-'},'patchSaturation',0.3);
        ylabel('Amplitude (uV)');
        
        % Axes and parameters
        title([tfr_avg{cond_ix}.label{ch_ix} ': ' cond_names{cond_ix}]);
        set(gca,'YLim', erp_ylim);
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        set(gca,'FontSize',16);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end

end
