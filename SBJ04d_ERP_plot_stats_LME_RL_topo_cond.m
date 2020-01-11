function SBJ04d_ERP_plot_stats_LME_RL_topo_cond(SBJs,proc_id,an_id,stat_id,plt_id,save_fig,varargin)
% Plots group RL beta topographies with significance for ERPs
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
[reg_lab, reg_colors, reg_styles]  = fn_regressor_label_styles(st.model_lab);
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles(st.trial_cond{1});

% Check for window compatibility
if strcmp(st.measure,'ts')
    error('Single topo plot must have single metric within window!');
end

%% Load Stats
load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' stat_id '_' an_id '.mat']);

% Get color limits
cfgat = [];
cfgat.latency = plt.plt_lim;
cfgat.avgovertime = 'yes';
clim = [0 0];
for cond_ix = 1:numel(cond_lab)
    tmp = ft_selectdata(cfgat,er_grp{cond_ix});
    clim = [min([clim(1) min(tmp.avg)]) max([clim(2) max(tmp.avg)])];
end


%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' stat_id '/' an_id '/' plt_id '/'];
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
    
    % Obtain Model Parameters
    r2 = NaN(size(st_time_vec));
    plot_betas = NaN([numel(reg_lab) numel(st_time_vec)]);
    for t_ix = 1:numel(st_time_vec)
        plot_betas(:,t_ix) = lme{t_ix}.Coefficients.Estimate(2:end);
        r2(t_ix) = lme{t_ix}.Rsquared.Adjusted;
    end
    
    % Find significant time periods
    sig_reg    = false(size(reg_lab));
    sig_chunks = cell(size(reg_lab));
    for reg_ix = 1:numel(reg_lab)
        if any(qvals(reg_ix,:) <= st.alpha)
            sig_reg(reg_ix) = true;
            if strcmp(st.measure,'ts')
                sig_chunks{reg_ix} = fn_find_chunks(squeeze(qvals(reg_ix,:))<=st.alpha);
                sig_chunks{reg_ix}(squeeze(qvals(reg_ix,sig_chunks{reg_ix}(:,1))>st.alpha),:) = [];
            elseif strcmp(st.measure,'mean')
                win_lim = zeros([1 2]);
                [~, win_lim(1)] = min(abs(time_vec - st.stat_lim(1)));
                [~, win_lim(2)] = min(abs(time_vec - st.stat_lim(2)));
                sig_chunks{reg_ix} = win_lim;
            end
        end
    end
    
    %% Create plot
    fig_name = ['GRP_' stat_id '_' ch_list{ch_ix}];
    if plot_median; fig_name = [fig_name '_med']; end
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 1],'Visible',fig_vis);   %this size is for single plots
    
    %% Plot ERPs with significance
    axes = gobjects([3 1]);
    subplot(3,1,1);
    axes(1) = gca; hold on;
    
    % Plot Means (and variance)
    cond_lines = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+sum(sig_reg)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        cond_lines{cond_ix} = shadedErrorBar(time_vec, plot_means(cond_ix,:), sems(cond_ix,:),...
            {'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},plt.errbar_alpha);
        main_lines(cond_ix) = cond_lines{cond_ix}.mainLine;
    end
    ylims = ylim;
    if strcmp(plt.sig_type,'line')
        data_lim = [min(min(plot_means-sems)) max(max(plot_means+sems))];
    end
    
    % Plot Significance
    for reg_ix = 1:numel(reg_lab)
        for sig_ix = 1:size(sig_chunks{reg_ix},1)
            if strcmp(plt.sig_type,'line')
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
                    main_lines(numel(cond_lab)+reg_ix) = sig_line;
                end
            elseif strcmp(plt.sig_type,'patch')
                if reg_ix>1; warning('Why use patch sig with more than 1 group???'); end
                patch(w2.time([sig_chunks{reg_ix}(sig_ix,1) sig_chunks{reg_ix}(sig_ix,1) ...
                    sig_chunks{reg_ix}(sig_ix,2) sig_chunks{reg_ix}(sig_ix,2)]),...
                    [ylims(1) ylims(2) ylims(2) ylims(1)],...
                    plt.sig_color,'FaceAlpha',plt.sig_alpha);
            elseif strcmp(plt.sig_type,'bold')
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
    %	Assume only one event at 0; ... evnt_ix = 1:numel(an.event_type)
    main_lines(end) = line([0 0],ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color{1},'LineStyle',plt.evnt_style{1});
    
    if strcmp(plt.sig_type,'line')
        leg_lab = [cond_lab reg_lab(sig_reg) an.event_type];
    else
        leg_lab = [cond_lab an.event_type];
    end
    
    % Axes and Labels
    axes(1).YLabel.String = 'uV';
    axes(1).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(1).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(1).XLabel.String = 'Time (s)';
    axes(1).Title.String  = ch_list{ch_ix};
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    ylims = ylim;
    set(gca,'FontSize',16);
    axes(1).YLim = ylims;
    
    %% Plot Betas and R2
    subplot(3,1,2);
    axes(2) = gca; hold on;
    
    % Plot Model Betas
    beta_lines = gobjects([numel(reg_lab)+numel(an.event_type) 1]);
    for reg_ix = 1:numel(reg_lab)
        beta_lines(reg_ix) = line(st_time_vec, plot_betas(reg_ix,:),...
            'Color',reg_colors{reg_ix},'LineWidth',plt.mean_width,...
            'LineStyle',reg_styles{reg_ix});
    
        % Plot Significance
        for sig_ix = 1:size(sig_chunks{reg_ix},1)
            % Assume: strcmp(plt.sig_type,'bold')
            sig_times = st_time_vec(sig_chunks{reg_ix}(sig_ix,1):sig_chunks{reg_ix}(sig_ix,2));
            line(sig_times,plot_betas(reg_ix,sig_chunks{reg_ix}(sig_ix,1):sig_chunks{reg_ix}(sig_ix,2)),...
                'Color',reg_colors{reg_ix},'LineStyle',reg_styles{reg_ix},...
                'LineWidth',plt.sig_width);
        end
    end
    
    % Plot Events
    %	Assume only one event at 0; ... evnt_ix = 1:numel(an.event_type)
    beta_lines(numel(reg_lab)+1) = line([0 0],ylim,...
        'LineWidth',plt.evnt_width,'Color',plt.evnt_color{1},'LineStyle',plt.evnt_style{1});
        
    % Axes and Labels
%     ax.YLim          = ylims; %!!! change for plt.sigType=line
    axes(2).YLabel.String = 'Beta Weight';
    axes(2).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(2).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(2).XLabel.String = 'Time (s)';
    axes(2).Title.String  = 'Model Weights';
    if plt.legend
        legend(beta_lines,[reg_lab an.event_type],'Location',plt.legend_loc);
    end
    ylims = ylim;
    set(gca,'FontSize',16);
    axes(2).YLim = ylims;
    
    %% Plot R2
    subplot(3,1,3);
    axes(3) = gca; hold on;
    
    line(st_time_vec, r2, 'Color','k', 'LineWidth',2);
    line([0 0],ylim, 'LineWidth',plt.evnt_width, 'Color',plt.evnt_color{1},...
        'LineStyle',plt.evnt_style{1});
    ylims = ylim;
    
    % Axes and Labels
    axes(3).YLabel.String = 'Adjusted R2';
    axes(3).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(3).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(3).XLabel.String = 'Time (s)';
    axes(3).Title.String  = 'Model Fit';
    set(gca,'FontSize',16);
    axes(3).YLim = ylims;
    
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
