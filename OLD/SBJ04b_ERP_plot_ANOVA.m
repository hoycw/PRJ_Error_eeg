function SBJ04b_ERP_plot_stats_ANOVA(SBJ,proc_id,an_id,plt_id,save_fig,varargin)
error('split into ANOVA and GLM!');
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
if ~exist('write_report','var'); write_report = 0; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load data
load([SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',an_id,'.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

% Select conditions (and trials)
[grp_lab, grp_colors, ~] = fn_group_label_styles(an.model_lab);
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles(an.model_lab);
full_cond_idx = fn_condition_index(cond_lab, bhv);
bhv_fields = fieldnames(bhv);
orig_n_trials = numel(bhv.trl_n);
for f_ix = 1:numel(bhv_fields)
    if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
        bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
    end
end
cond_idx = fn_condition_index(cond_lab, bhv);

% Get trials for plotting
trials = cell(size(cond_lab));
cond_idx = fn_condition_index(cond_lab, bhv);
for cond_ix = 1:numel(cond_lab)
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(roi.label) numel(cond_trial_ix) numel(roi.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
    end
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' SBJ '/' an_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = zeros([numel(w2.label) numel(an.groups)]);
for ch_ix = 1:numel(w2.label)
    %% Compute plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(roi.time{1})]);
    sems  = NaN([numel(cond_lab) numel(roi.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(ch_ix,:,:),2));
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(ch_ix,:,:),[],2))./sqrt(size(trials{cond_ix},2))';
    end
    
    % Find significant time periods
    sig_grp = false(size(an.groups));
    sig_chunks = cell(size(an.groups));
    for grp_ix = 1:numel(an.groups)
        if any(w2.qval(grp_ix,ch_ix,:) <= an.alpha)
            sig_grp(grp_ix) = true;
            sig_ch(ch_ix,grp_ix) = 1;
            sig_chunks{grp_ix} = fn_find_chunks(squeeze(w2.qval(grp_ix,ch_ix,:))<=an.alpha);
            sig_chunks{grp_ix}(squeeze(w2.qval(grp_ix,ch_ix,sig_chunks{grp_ix}(:,1))>an.alpha),:) = [];
        end
    end
    
    %% Create plot
    fig_name = [SBJ '_' an_id '_' w2.label{ch_ix}];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
%     [plot_rc,~] = fn_num_subplots(numel(w2.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    ax = gca; hold on;
    
    % Plot individual trials per condition
    if plt.butterfly
        for cond_ix = 1:numel(cond_lab)
            plot(roi.time{1},squeeze(trials{cond_ix}(ch_ix,:,:)),...
                'Color',cond_colors{cond_ix},'LineWidth',plt.butterfly_width,...
                'LineStyle',cond_styles{cond_ix});
        end
    end

    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(roi.time{1}, means(cond_ix,:), sems(cond_ix,:),...
            {'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    ylims = ylim;
    if strcmp(plt.sig_type,'line')
        data_lim = [min(min(means-sems)) max(max(means+sems))];
    end
    
    % Plot Extra Features (events, significance)
    for evnt_ix = 1:numel(an.event_type)
        main_lines(numel(cond_lab)+evnt_ix) = line([0 0],ylim,...
            'LineWidth',plt.evnt_width(evnt_ix),'Color',plt.evnt_color{evnt_ix},'LineStyle',plt.evnt_style{evnt_ix});
    end
    
    % Plot Significance
    for grp_ix = 1:numel(an.groups)
        for sig_ix = 1:size(sig_chunks{grp_ix},1)
            if strcmp(plt.sig_type,'line')
                sig_times = w2.time(sig_chunks{grp_ix}(sig_ix,1):sig_chunks{grp_ix}(sig_ix,2));
                if strcmp(plt.sig_loc,'below')
                    sig_y = ylims(1) - grp_ix*data_lim(1)*plt.sig_loc_factor;
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
    ax.Title.String  = w2.label{ch_ix};
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_filename]);
    end
end

%% Save out list of channels with significant differences
if write_report
    sig_report_fname = [fig_dir 'ch_sig_list.csv'];
    sig_report = fopen(sig_report_fname,'w');
    fprintf(sig_report,'%s\n',an_id);
    fprintf(sig_report,'label,%s\n',an.model_lab);
    for ch_ix = 1:numel(w2.label)
        fprintf(sig_report,'%s,%.0f\n',w2.label{ch_ix},sig_ch(ch_ix));
    end
    fclose(sig_report);
end

end
