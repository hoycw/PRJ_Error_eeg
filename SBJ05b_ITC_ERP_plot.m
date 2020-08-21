function SBJ05b_ITC_ERP_plot(SBJ,conditions,proc_id,itc_an_id,erp_an_id,plt_id,save_fig,varargin)
%% Compute and plot ITPC matrix per condition for single SBJ; plot ERP on top
% INPUTS:
%   SBJ [str] - ID of subject to plot
%   conditions [str] - group of condition labels to segregate trials
%   proc_id [str] - ID of preprocessing pipeline
%   itc_an_id [str] - ID of the TFR/ITPC analysis parameters to use
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
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' erp_an_id '_vars.m'];
eval(an_vars_cmd);
erp_an = an;
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' itc_an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~an.complex; error('why run this without ITPC an_vars?'); end
if ~strcmp(an.event_type,erp_an.event_type); error('itc and erp event mismatch!'); end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load data
load([SBJ_vars.dirs.proc SBJ '_' erp_an_id '.mat']);
load([SBJ_vars.dirs.proc SBJ '_' proc_id '_' itc_an_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

% Select conditions (and trials)
[grp_lab, ~, ~, ~] = fn_group_label_styles(conditions);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(conditions);

% Group conditions for subplot organization
grp_cond_lab = cell(size(grp_lab));
for grp_ix = 1:numel(grp_lab)
    [grp_cond_lab{grp_ix}, ~, ~, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
end
cond_idx = fn_condition_index(cond_lab, bhv);

% Get trials for ERP plotting
trials = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(roi.label) numel(cond_trial_ix) numel(roi.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
    end
end

%% Compute Inter-Trial Phase Coherence
itpc = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    % Compute ITPC
    F = tfr.fourierspctrm(cond_idx==cond_ix,:,:,:);
    itpc{cond_ix} = F./abs(F);                                  % Normalize to unit circle
    itpc{cond_ix} = sum(itpc{cond_ix},1);                       % Sum phase angles
    itpc{cond_ix} = abs(itpc{cond_ix})/sum(cond_idx==cond_ix);	% Get mean of angles for consistency
end

%% Get event timing for plotting
if strcmp(an.event_type,'S')
    error('add loading of prdm_vars to plot relative to stim!');
end
[evnt_times] = fn_get_evnt_times(an.event_type,plt.evnt_lab);

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' itc_an_id '/' conditions '/' erp_an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr.label)
    %% Compute ERP plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(roi.time{1})]);
    sems  = NaN([numel(cond_lab) numel(roi.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(ch_ix,:,:),2));
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(ch_ix,:,:),[],2))./sqrt(size(trials{cond_ix},2))';
    end
    
    %% Create plot
    fig_name = [SBJ '_' conditions '_' itc_an_id '_' erp_an_id '_' tfr.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Get color lims per condition
    clim = zeros([numel(cond_lab) 2]);
    for cond_ix = 1:numel(cond_lab)
        vals = itpc{cond_ix}(1,ch_ix,:,:);
        clim(cond_ix,:) = [min(vals(:)) max(vals(:))];
    end
    tick_ix = 1:3:numel(tfr.freq);
    yticklab = cell(size(tick_ix));
    for f = 1:numel(tick_ix)
        yticklab{f} = num2str(tfr.freq(tick_ix(f)),'%.1f');
    end
    
    % Condition Plots
    for cond_ix = 1:length(cond_lab)
        subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix);
        
        % Plot ITC Matrix
        yyaxis left
        %contourf(tfr.time, tfr.freq, squeeze(itpc{cond_ix}(1,ch_ix,:,:)));
        imagesc(tfr.time, tfr.freq, squeeze(itpc{cond_ix}(1,ch_ix,:,:)));% 1:numel(tfr.freq)
        set(gca,'YDir','normal');
%         set(gca,'YTick',1:3:numel(tfr.freq));
%         set(gca,'YTickLabels',yticklab);
        ylabel('Frequency (Hz)');
        caxis([min(clim(:,1)) max(clim(:,2))]);
        colorbar('northoutside');
        
        % Plot Events
        for evnt_ix = 1:numel(plt.evnt_lab)
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        % Plot ERP Means (and variance)
        yyaxis right
        shadedErrorBar(roi.time{1}, means(cond_ix,:), sems(cond_ix,:),...
                'lineProps',{'Color','k','LineWidth',2,...
                'LineStyle','-'},'patchSaturation',0.3);
        ylabel('Amplitude (uV)');
        
        % Axes and parameters
        title([tfr.label{ch_ix} ': ' cond_lab{cond_ix}]);
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        set(gca,'FontSize',16);
    end
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
