function SBJ06f_CPA_candidate_ERP_plot_GRP_betas(SBJs,eeg_proc_id,cpa_id,an_id,stat_id,plt_id,save_fig,varargin)
% Plots group CPA candidate betas after SBJ stats, one per regressor
%   Only for single channel right now...
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

% Select Conditions of Interest
[reg_lab, ~, reg_colors, reg_styles]  = fn_regressor_label_styles(st.model_lab);

%% Load ERPs
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %load([SBJ_vars.dirs.events SBJs{s} '_behav_' eeg_proc_id '_final.mat']);
    
    % Get time and channel info
    if s==1
        load([SBJ_vars.dirs.proc SBJs{s} '_' cpa_id '_' an_id '.mat']);
        time_vec = roi.time{1};
        ch_list  = roi.label;
        if numel(ch_list)>1; error('More than 1 channel, not ready!'); end
        
        % Select time and trials of interest
        cfgs = [];
        cfgs.latency = st.stat_lim;
        st_roi = ft_selectdata(cfgs, roi);
        st_time_vec = st_roi.time{1};
        
        betas = nan([numel(reg_lab) numel(SBJs) numel(ch_list) numel(st_time_vec)]);
        r2    = nan([numel(SBJs) numel(ch_list) numel(st_time_vec)]);
    end
    
    % Load beta and R2 values
    load([SBJ_vars.dirs.proc SBJs{s} '_' stat_id '_' cpa_id '_' an_id '.mat']);
    for t_ix = 1:numel(st_time_vec)
        betas(:,s,1,t_ix) = lme{t_ix}.Coefficients.Estimate(2:end);
        r2(s,1,t_ix) = lme{t_ix}.Rsquared.Adjusted;
    end
    
    clear lme qvals SBJ_vars roi st_roi
end

% %% Group stats
% pvals = NaN([numel(reg_lab) numel(st_time_vec)]);
% for reg_ix = 1:numel(reg_lab)
%     
% end

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
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/candidate/' cpa_id '_' an_id '/' stat_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Compute means and variance
    plot_betas = NaN([numel(reg_lab) numel(st_time_vec)]);
    beta_sems  = NaN([numel(reg_lab) numel(st_time_vec)]);
    for reg_ix = 1:numel(reg_lab)
        if plot_median
            plot_betas(reg_ix,:) = squeeze(median(betas(reg_ix,:,ch_ix,:),2));
        else
            plot_betas(reg_ix,:) = squeeze(mean(betas(reg_ix,:,ch_ix,:),2));
        end
        beta_sems(reg_ix,:) = squeeze(std(betas(reg_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    if plot_median
        plot_r2 = squeeze(median(r2(:,ch_ix,:),1));
    else
        plot_r2 = squeeze(mean(r2(:,ch_ix,:),1));
    end
    r2_sem = squeeze(std(r2(:,ch_ix,:),[],1))./sqrt(numel(SBJs))';
    
    % Get beta limits
    ch_betas = squeeze(betas(:,:,ch_ix,:));
    beta_lim = [min(ch_betas(:)) max(ch_betas(:))];
    
    %% Create plot
    fig_name = ['GRP_betas_' stat_id '_' cpa_id '_' an_id '_' ch_list{ch_ix}];
    if plot_median; fig_name = [fig_name '_med']; end
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);   %this size is for single plots
    
    %% Plot betas
    axes = gobjects([numel(reg_lab)+1 1]);
    [num_rc,~] = fn_num_subplots(numel(reg_lab)+1);
    for reg_ix = 1:numel(reg_lab)
        subplot(num_rc(1),num_rc(2),reg_ix);
        axes(reg_ix) = gca; hold on;
        
        % Plot individual betas
        if plt.butterfly
            plot(st_time_vec,squeeze(betas(reg_ix,:,ch_ix,:)),...
                'Color','k','LineWidth',plt.butterfly_width,...
                'LineStyle',reg_styles{reg_ix});
        end
        
        % Plot mean betas
        main_lines = gobjects([numel(plt.evnt_lab)+1 1]);
        r2_line = shadedErrorBar(st_time_vec, plot_betas(reg_ix,:), beta_sems(reg_ix,:),...
            {'Color',reg_colors{reg_ix},'LineWidth',plt.mean_width,...
            'LineStyle',reg_styles{reg_ix}},plt.errbar_alpha);
        main_lines(1) = r2_line.mainLine;
        
        % Plot Events
        for evnt_ix = 1:numel(plt.evnt_lab)
            main_lines(evnt_ix+1) = line(...
                [evnt_times(evnt_ix) evnt_times(evnt_ix)],beta_lim,...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        leg_lab = [{'Mean beta'} plt.evnt_lab];
        
        % Axes and Labels
        axes(reg_ix).YLabel.String = 'Beta';
        axes(reg_ix).YLim          = beta_lim;
        axes(reg_ix).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
        axes(reg_ix).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
        axes(reg_ix).XLabel.String = 'Time (s)';
        axes(reg_ix).Title.String  = [ch_list{ch_ix} ': ' reg_lab{reg_ix}];
        if plt.legend
            legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
        end
        ylims = ylim;
        set(gca,'FontSize',16);
        axes(reg_ix).YLim = ylims;
    end
    
    %% Plot R2
    subplot(num_rc(1),num_rc(2),numel(reg_lab)+1);
    axes(end) = gca; hold on;
    
    % Plot individual R2s
    if plt.butterfly
        plot(st_time_vec,squeeze(r2(:,ch_ix,:)),...
            'Color','k','LineWidth',plt.butterfly_width,...
            'LineStyle','-');
    end
    
    % Plot mean R2
    main_lines = gobjects([numel(plt.evnt_lab)+1 1]);
    r2_line = shadedErrorBar(st_time_vec, plot_r2, r2_sem,...
            {'Color',reg_colors{reg_ix},'LineWidth',plt.mean_width,...
            'LineStyle',reg_styles{reg_ix}},plt.errbar_alpha);
        main_lines(1) = r2_line.mainLine;
    main_lines(1) = line(st_time_vec, squeeze(mean(r2(:,ch_ix,:),1)), 'Color','k',...
        'LineWidth',plt.mean_width);
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_lines(evnt_ix+1) = line(...
            [evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end
    
    leg_lab = [{'Mean beta'} plt.evnt_lab];
    
    % Axes and Labels
    axes(end).YLabel.String = 'R2';
    axes(end).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(end).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(end).XLabel.String = 'Time (s)';
    axes(end).Title.String  = [ch_list{ch_ix} ': Model Fit'];
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    ylims = ylim;
    set(gca,'FontSize',16);
    axes(end).YLim = ylims;
    
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
