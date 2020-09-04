function SBJ04d_ERP_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_id,stat_id,plt_id,save_fig,varargin)
%% Plots group ERPs with significance, model coefficients, and model fit
%   Grand average ERPs with significant epochs marked underneath
%       Optional: grand median instead of grand average
%   Beta weight time series per regressor, bolded for significant epochs
%   R2 time series for model fit
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
%   plt_id [str] - ID of the plotting parameters to use
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
%       plot_median [0/1] - binary flag to plot median ERP instead of mean
%           default: 0
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
        elseif strcmp(varargin{v},'plot_median')
            plot_median = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('plot_median','var'); plot_median = 0; end
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

%% Load Stats
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat']);

% Load paradigm variables to plot event timing
warning('WARNING: Assuming same prdm_vars for all SBJ to get event timing!');
prdm_vars = load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/03_events/' SBJs{1} '_prdm_vars.mat']);

%% Load ERPs
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    load([SBJ_vars.dirs.proc,SBJs{s},'_',an_id,'.mat']);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    
    if s==1
        time_vec = roi.time{1};
        ch_list  = roi.label;
        means = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
        
        % Select time and trials of interest
        cfgs = [];
        cfgs.latency = st.stat_lim;
        st_roi = ft_selectdata(cfgs, roi);
        st_time_vec = st_roi.time{1};
    end
    
    % Select Conditions of Interest
    cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for t_ix = 1:numel(cond_trial_ix)
            trials(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
        end
        
        % Average ERP within condition per SBJ
        means(cond_ix,s,:,:) = mean(trials,2);
    end
    
    clear tmp SBJ_vars bhv roi
end

%% Get event timing for plotting
[evnt_times] = fn_get_evnt_times(an.event_type,plt.evnt_lab,'prdm_vars',prdm_vars);

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' stat_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Compute means and variance
    plot_means = NaN([numel(cond_lab) numel(time_vec)]);
    sems  = NaN([numel(cond_lab) numel(time_vec)]);
    for cond_ix = 1:numel(cond_lab)
        if plot_median
            plot_means(cond_ix,:) = squeeze(median(means(cond_ix,:,ch_ix,:),2));
        else
            plot_means(cond_ix,:) = squeeze(mean(means(cond_ix,:,ch_ix,:),2));
        end
        sems(cond_ix,:) = squeeze(std(means(cond_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    
    % Obtain Model Coefficients and Fit
    r2 = NaN(size(st_time_vec));
    plot_betas = NaN([numel(reg_lab) numel(st_time_vec)]);
    for t_ix = 1:numel(st_time_vec)
        if strcmp(st.model_lab,'SBJonly')
            plot_betas(:,t_ix) = lme{t_ix}.Coefficients.Estimate;
        else
            plot_betas(:,t_ix) = lme{t_ix}.Coefficients.Estimate(2:end);
        end
        r2(t_ix) = lme{t_ix}.Rsquared.Adjusted;
    end
    
    % Find significant time periods
    sig_reg    = false(size(reg_lab));
    sig_chunks = cell(size(reg_lab));
    for reg_ix = 1:numel(reg_lab)
        if any(qvals(reg_ix,:) <= st.alpha)
            sig_reg(reg_ix) = true;
            if strcmp(st.measure,'ts')
                % Find consecutive chunks of (non-)significant coefficients
                sig_chunks{reg_ix} = fn_find_chunks(squeeze(qvals(reg_ix,:))<=st.alpha);
                % Remove non-significant chunks
                sig_chunks{reg_ix}(squeeze(qvals(reg_ix,sig_chunks{reg_ix}(:,1))>st.alpha),:) = [];
            elseif strcmp(st.measure,'mean')
                % Find window edges
                win_lim = zeros([1 2]);
                [~, win_lim(1)] = min(abs(time_vec - st.stat_lim(1)));
                [~, win_lim(2)] = min(abs(time_vec - st.stat_lim(2)));
                sig_chunks{reg_ix} = win_lim;
            end
        end
    end
    
    %% Create plot
    fig_name = [SBJ_id '_' stat_id '_' ch_list{ch_ix}];
    if plot_median; fig_name = [fig_name '_med']; end
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.7 1],'Visible',fig_vis);
    
    %% Plot ERPs with significance
    axes = gobjects([3 1]);
    subplot(6,1,1:3);
    axes(1) = gca; hold on;
    
    % Plot ERP Means (and variance)
    cond_lines = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+sum(sig_reg)+numel(plt.evnt_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        cond_lines{cond_ix} = shadedErrorBar(time_vec, plot_means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(cond_ix) = cond_lines{cond_ix}.mainLine;
    end
    
    % Fix y limits to be consistent across channels
    if any(strcmp(SBJ_id,{'good1','goodall'}))
        ylims = [-15 30];
    else
        ylims = ylim;
    end
    
    % Find min/max error bars to plot significant epochs underneath
    if strcmp(plt.sig_type,'line')
        data_lim = [min(min(plot_means-sems)) max(max(plot_means+sems))];
    end
    
    % Plot Significance per Regressor
    sig_reg_ix = find(sig_reg);
    for r = 1:numel(sig_reg_ix)
        reg_ix = sig_reg_ix(r);
        for sig_ix = 1:size(sig_chunks{reg_ix},1)
            if strcmp(plt.sig_type,'line')
                % Plot horizontal line beneath/above the ERPs
                sig_times = st_time_vec(sig_chunks{reg_ix}(sig_ix,1):sig_chunks{reg_ix}(sig_ix,2));
                if strcmp(plt.sig_loc,'below')
                    sig_y = data_lim(1) + reg_ix*data_lim(1)*plt.sig_loc_factor;
                elseif strcmp(plt.sig_loc,'above')
                    sig_y = data_lim(2) + reg_ix*data_lim(2)*plt.sig_loc_factor;
                end
                sig_line = line(sig_times,repmat(sig_y,size(sig_times)),...
                    'LineWidth',plt.sig_width,'Color',reg_colors{reg_ix},...
                    'LineStyle',reg_styles{reg_ix});
                if sig_ix==1
                    main_lines(numel(cond_lab)+r) = sig_line;
                end
            elseif strcmp(plt.sig_type,'patch')
                % Plot rectangular patch over epoch
                if reg_ix>1; warning('Why use patch sig with more than 1 group???'); end
                patch(w2.time([sig_chunks{reg_ix}(sig_ix,1) sig_chunks{reg_ix}(sig_ix,1) ...
                    sig_chunks{reg_ix}(sig_ix,2) sig_chunks{reg_ix}(sig_ix,2)]),...
                    [ylims(1) ylims(2) ylims(2) ylims(1)],...
                    plt.sig_color,'FaceAlpha',plt.sig_alpha);
            elseif strcmp(plt.sig_type,'bold')
                % Bold the ERP time series
                error('bold sig type not ready');
%                 sig_times = sig_chunks{grp_ix}(sig_ix,1):sig_chunks{grp_ix}(sig_ix,2);
%                 for cond_ix = 1:numel(cond_lab)
%                     line(sig_times,means(ch_ix,sig_chunks{grp_ix}(sig_ix,1):sig_chunks{grp_ix}(sig_ix,2)),...
%                         'Color',cond_colors{cond_ix},'LineStyle',plt.sig_style,...
%                         'LineWidth',plt.sig_width);
%                 end
            else
                error('unknown sig_type');
            end
        end
    end
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_lines(numel(cond_lab)+sum(sig_reg)+evnt_ix) = line(...
            [evnt_times(evnt_ix) evnt_times(evnt_ix)],[-15 30],...%ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end    
    if strcmp(plt.sig_type,'line')
        leg_lab = [cond_lab reg_lab(sig_reg) plt.evnt_lab];
    else
        leg_lab = [cond_lab plt.evnt_lab];
    end
    
    % Axes and Labels
    axes(1).YLabel.String = 'uV';
    axes(1).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(1).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(1).XLabel.String = 'Time (s)';
    axes(1).Title.String  = [ch_list{ch_ix} ' (n = ' num2str(numel(SBJs)) ')'];
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    % Fix y limits to be consistent across channels
    if any(strcmp(SBJ_id,{'good1','goodall'}))
        ylims = [-15 30];
    else
        ylims = ylim;
    end
    set(gca,'FontSize',16);
    axes(1).YLim = ylims;
    
    %% Plot Betas
    subplot(6,1,4:5);
    axes(2) = gca; hold on;
    
    % Plot Model Betas
    beta_lines = gobjects([numel(reg_lab)+numel(plt.evnt_lab) 1]);
    for reg_ix = 1:numel(reg_lab)
        beta_lines(reg_ix) = line(st_time_vec, plot_betas(reg_ix,:),...
            'Color',reg_colors{reg_ix},'LineWidth',plt.mean_width,...
            'LineStyle',reg_styles{reg_ix});
    
        % Plot Significance by bolding time series
        for sig_ix = 1:size(sig_chunks{reg_ix},1)
            % Assume: strcmp(plt.sig_type,'bold')
            sig_times = st_time_vec(sig_chunks{reg_ix}(sig_ix,1):sig_chunks{reg_ix}(sig_ix,2));
            line(sig_times,plot_betas(reg_ix,sig_chunks{reg_ix}(sig_ix,1):sig_chunks{reg_ix}(sig_ix,2)),...
                'Color',reg_colors{reg_ix},'LineStyle',reg_styles{reg_ix},...
                'LineWidth',plt.sig_width);
        end
    end
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        beta_lines(numel(reg_lab)+evnt_ix) = line(...
            [evnt_times(evnt_ix) evnt_times(evnt_ix)],[-4 6],...%ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end    
        
    % Axes and Labels
%     ax.YLim          = ylims; %!!! change for plt.sigType=line
    axes(2).YLabel.String = 'Beta Weight';
    axes(2).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(2).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(2).XLabel.String = 'Time (s)';
    % axes(2).Title.String  = 'Model Weights';
    if plt.legend
        legend(beta_lines,[reg_lab plt.evnt_lab],'Location',plt.legend_loc);
    end
    % Fix y limits to be consistent across channels
    if strcmp(SBJ_id,'goodall')
        ylims = [-3 5];
    elseif strcmp(SBJ_id,'good1')
        ylims = [-4 6];
    else
        ylims = ylim;
    end
    set(gca,'FontSize',16);
    axes(2).YLim = ylims;
    
    %% Plot Model Fit
    subplot(6,1,6);
    axes(3) = gca; hold on;
    
    % Plot R2
    line(st_time_vec, r2, 'Color','k', 'LineWidth',2);
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        line([evnt_times(evnt_ix) evnt_times(evnt_ix)],[0 0.4],...%ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end    
    ylims = ylim;
    
    % Axes and Labels
    axes(3).YLabel.String = 'Adjusted R2';
    axes(3).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(3).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(3).XLabel.String = 'Time (s)';
    % axes(3).Title.String  = 'Model Fit';
    set(gca,'FontSize',16);
    axes(3).YLim = ylims;
    
    %% Report peak stats per regressor
    % Prints largest model coefficient, time point, and q value
    for reg_ix = 1:numel(reg_lab)
        max_tmp = max(plot_betas(reg_ix,:));
        min_tmp = min(plot_betas(reg_ix,:));
        if abs(max_tmp) > abs(min_tmp)
            [max_beta, max_t_ix] = max(plot_betas(reg_ix,:));
        else
            [max_beta, max_t_ix] = min(plot_betas(reg_ix,:));
        end
        fprintf('%s max beta = %.03f at %.03f; p = %.20f\n',reg_lab{reg_ix},max_beta,...
            st_time_vec(max_t_ix),qvals(reg_ix,max_t_ix));
    end
    
    % Print maximum model fit, time point
    [max_r2, max_t_ix] = max(r2);
    fprintf('max R2 = %.03f at %.03f\n',max_r2,st_time_vec(max_t_ix));
    
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

end
