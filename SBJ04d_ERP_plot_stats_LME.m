function SBJ04d_ERP_plot_stats_LME(SBJs,proc_id,an_id,stat_id,plt_id,save_fig,varargin)
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
[grp_lab, grp_colors, ~]  = fn_group_label_styles(st.model_lab);
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles(st.model_lab);
fact_cond_lab = cell(size(st.factors));
for fct_ix = 1:numel(st.factors)
    [fact_cond_lab{fct_ix}, ~, ~, ~] = fn_condition_label_styles(st.factors{fct_ix});
    if numel(fact_cond_lab{fct_ix}) > 2; error('coding not ready for 3 conditions'); end
end

%% Load Stats
load([root_dir 'PRJ_Error_eeg/data/GRP/' 'GRP_' stat_id '_' an_id '.mat']);

% Adjust p values
pvals = nan([numel(grp_lab) numel(lme)]);
qvals = nan([numel(grp_lab) numel(lme)]);
for grp_ix = 1:numel(grp_lab)
    label_ix = strcmp(lme{1}.CoefficientNames,strrep(grp_lab{grp_ix},'*',':'));
    for t_ix = 1:numel(lme)
        pvals(grp_ix,t_ix) = lme{t_ix}.Coefficients.pValue(label_ix);
    end
end

[~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);

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
    end
    
    % Select Conditions of Interest
    cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for t_ix = 1:numel(cond_trial_ix)
            trials(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
        end
        means(cond_ix,s,:,:) = mean(trials,2);
    end
    
    clear tmp SBJ_vars bhv roi
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
    
    % Find significant time periods
    sig_grp    = false(size(grp_lab));
    sig_chunks = cell(size(grp_lab));
    for grp_ix = 1:numel(grp_lab)
        if any(qvals(grp_ix,:) <= st.alpha)
            sig_grp(grp_ix) = true;
            sig_chunks{grp_ix} = fn_find_chunks(squeeze(qvals(grp_ix,:))<=st.alpha);
            sig_chunks{grp_ix}(squeeze(qvals(grp_ix,sig_chunks{grp_ix}(:,1))>st.alpha),:) = [];
        end
    end
    
    %% Create plot
    fig_name = ['GRP_' stat_id '_' ch_list{ch_ix}];
    if plot_median; fig_name = [fig_name '_med']; end
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
%     [plot_rc,~] = fn_num_subplots(numel(w2.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    ax = gca; hold on;
    
    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(time_vec, plot_means(cond_ix,:), sems(cond_ix,:),...
            {'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    ylims = ylim;
    if strcmp(plt.sig_type,'line')
        data_lim = [min(min(plot_means-sems)) max(max(plot_means+sems))];
    end
    
    % Plot Significance
    for grp_ix = 1:numel(grp_lab)
        for sig_ix = 1:size(sig_chunks{grp_ix},1)
            if strcmp(plt.sig_type,'line')
                [~, stat_start] = min(abs(time_vec-st.stat_lim(1)));
                sig_times = time_vec([sig_chunks{grp_ix}(sig_ix,1):sig_chunks{grp_ix}(sig_ix,2)]+stat_start);
                if strcmp(plt.sig_loc,'below')
                    sig_y = ylims(1) + grp_ix*data_lim(1)*plt.sig_loc_factor;
                elseif strcmp(plt.sig_loc,'above')
                    sig_y = ylims(2) + grp_ix*data_lim(2)*plt.sig_loc_factor;
                end
                sig_line = line(sig_times,repmat(sig_y,size(sig_times)),...
                    'LineWidth',plt.sig_width,'Color',grp_colors{grp_ix});
                if sig_ix==1
                    main_lines(end+1) = sig_line;
                end
            elseif strcmp(plt.sig_type,'patch')
                if grp_ix>1; warning('Why use patch sig with more than 1 group???'); end
                patch(w2.time([sig_chunks{grp_ix}(sig_ix,1) sig_chunks{grp_ix}(sig_ix,1) ...
                    sig_chunks{grp_ix}(sig_ix,2) sig_chunks{grp_ix}(sig_ix,2)]),...
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
    for evnt_ix = 1:numel(an.event_type)
        main_lines(numel(cond_lab)+evnt_ix) = line([0 0],ylim,...
            'LineWidth',plt.evnt_width(evnt_ix),'Color',plt.evnt_color{evnt_ix},'LineStyle',plt.evnt_style{evnt_ix});
    end
    
    if strcmp(plt.sig_type,'line')
        leg_lab = [cond_lab an.event_type grp_lab(sig_grp)];
    else
        leg_lab = [cond_lab an.event_type];
    end
    
    % Axes and Labels
%     ax.YLim          = ylims; %!!! change for plt.sigType=line
    ax.YLabel.String = 'uV';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    ax.Title.String  = ch_list{ch_ix};
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
