function SBJ05c_TFR_ERP_plot_grp(SBJs,conditions,proc_id,tfr_an_id,erp_an_id,plt_id,save_fig,varargin)
%% Plot TFRs averaged across the group
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
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' erp_an_id '_vars.m'];
eval(an_vars_cmd);
erp_an = an;
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' tfr_an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~strcmp(an.event_type,erp_an.event_type); error('itc and erp event mismatch!'); end
if ~all(an.trial_lim_s==erp_an.trial_lim_s); error('itc and erp trial_lim_s mismatch!'); end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select conditions (and trials)
[grp_lab, ~, ~, ~] = fn_group_label_styles(conditions);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(conditions);
% if ~strcmp(st.model_lab,{'DifOut','Out'}); error('not ready for surprise trials!'); end
grp_cond_lab = cell(size(grp_lab));
for grp_ix = 1:numel(grp_lab)
    [grp_cond_lab{grp_ix}, ~, ~, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
end

% Load data
tfr_all = cell([numel(cond_lab) numel(SBJs)]);
for s = 1:numel(SBJs)
    fprintf('-------------------- Processing %s ------------------------\n',SBJs{s});
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    load([SBJ_vars.dirs.proc SBJs{s} '_' erp_an_id '.mat'],'roi');
    load([SBJ_vars.dirs.proc SBJs{s} '_' proc_id '_' tfr_an_id '.mat'],'tfr');
    
    if s==1
        roi_time_vec = roi.time{1};
        means    = nan([numel(cond_lab) numel(SBJs) numel(tfr.label) numel(roi_time_vec)]);
    end
    
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
evnt_times = zeros(size(plt.evnt_lab));
if strcmp(an.event_type,'S')
    for evnt_ix = 1:numel(plt.evnt_lab)
        switch plt.evnt_lab{evnt_ix}
            case 'S'
                evnt_times(evnt_ix) = 0;
            case 'R'
                evnt_times(evnt_ix) = prdm_vars.target;
            case {'F','Fon'}
                evnt_times(evnt_ix) = prdm_vars.target+prdm_vars.fb_delay;
            case 'Foff'
                evnt_times(evnt_ix) = prdm_vars.target+prdm_vars.fb_delay+prdm_vars.fb;
            otherwise
                error(['Unknown event type in plt: ' plt.evnt_lab{evnt_ix}]);
        end
    end
elseif strcmp(an.event_type,'F')
    evnt_times(1) = 0;
else
    error('Unknown an.event_type');
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' tfr_an_id '/' conditions '/' erp_an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(tfr_avg{1}.label)
    %% Compute plotting data    
    % Compute means and variance
    plt_means = NaN([numel(cond_lab) numel(roi_time_vec)]);
    sems  = NaN([numel(cond_lab) numel(roi_time_vec)]);
    for cond_ix = 1:numel(cond_lab)
        plt_means(cond_ix,:) = squeeze(mean(means(cond_ix,:,ch_ix,:),2));
        sems(cond_ix,:) = squeeze(std(means(cond_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    
    %% Create plot
    fig_name = ['GRP_' conditions '_' tfr_an_id '_' erp_an_id '_' tfr_avg{1}.label{ch_ix}];
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
        yyaxis left
        imagesc(tfr_all{cond_ix}.time, tfr_all{cond_ix}.freq, squeeze(tfr_all{cond_ix}.powspctrm(ch_ix,:,:)));%1:numel(tfr_all{cond_ix}.freq)
        set(gca,'YDir','normal');
%         set(gca,'YTick',1:3:numel(tfr_all{cond_ix}.freq));
%         set(gca,'YTickLabels',yticklab);
        caxis([min(clim(:,1)) max(clim(:,2))]);
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
        
        title([tfr_avg{cond_ix}.label{ch_ix} ': ' cond_lab{cond_ix}]);
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
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
