function SBJ03c_ERP_plot_group(SBJs,proc_id,an_id,plt_id,save_fig,varargin)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

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
%         elseif strcmp(varargin{v},'write_report')
%             write_report = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
% if ~exist('write_report','var'); write_report = 0; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

[grp_lab, grp_colors, ~] = fn_group_label_styles(an.model_lab);
[cond_lab, cond_colors, cond_style, ~] = fn_condition_label_styles(an.model_lab);

% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m'];
eval(SBJ_vars_cmd);
tmp = load([SBJ_vars.dirs.SBJ,'04_proc/',SBJs{1},'_',an_id,'.mat'],'roi');
SBJs_vars = cell(size(SBJs));
rois  = cell(size(SBJs));
means = nan([numel(cond_lab) numel(SBJs) numel(roi.label) numel(tmp.roi.time{1})]);
w2s   = cell(size(SBJs));
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    SBJs_vars{s} = SBJ_vars;
    
    load([SBJs_vars{s}.dirs.proc,SBJs{s},'_',an_id,'.mat']);
    w2s{s} = w2; rois{s} = roi;
    load([SBJs_vars{s}.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    
    
    % Select Conditions of Interest
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    bhv_fields = fieldnames(bhv);
    orig_n_trials = numel(bhv.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
            bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
        end
    end
    
    cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        trials = nan([numel(roi.label) numel(cond_trial_ix) numel(roi.time{1})]);
        for t_ix = 1:numel(cond_trial_ix)
            trials(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
        end
        means(cond_ix,s,:,:) = mean(trials,2);
    end
    
    clear tmp SBJ_vars w2 roi
end

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
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/GRP/' an_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = zeros(size(w2{1}.label));
for ch_ix = 1:numel(w2.label)
    %% Compute plotting data    
    % Compute means and variance
    plt_means = NaN([numel(cond_lab) numel(roi.time{1})]);
    sems  = NaN([numel(cond_lab) numel(roi.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        plt_means(cond_ix,:) = squeeze(mean(means(cond_ix,:,ch_ix,:),2));
        sems(cond_ix,:) = squeeze(std(means(cond_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
    end
    
%     % Find significant time periods
%     sig_grp = false(size(an.groups));
%     sig_chunks = cell(size(an.groups));
%     for grp_ix = 1:numel(an.groups)
%         if any(w2.qval(grp_ix,ch_ix,:) <= an.alpha)
%             sig_grp(grp_ix) = true;
%             sig_ch(ch_ix,grp_ix) = 1;
%             sig_chunks{grp_ix} = fn_find_chunks(squeeze(w2.qval(grp_ix,ch_ix,:))<=an.alpha);
%             sig_chunks{grp_ix}(squeeze(w2.qval(grp_ix,ch_ix,sig_chunks{grp_ix}(:,1))>an.alpha),:) = [];
%         end
%     end
    
    %% Create plot
    fig_name = [SBJ '_' an_id '_' w2.label{ch_ix}];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
%     [plot_rc,~] = fn_num_subplots(numel(w2.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    ax = gca; hold on;
    
%     % Plot individual trials per condition
%     if plt.butterfly
%         for cond_ix = 1:numel(cond_lab)
%             plot(roi.time{1},squeeze(trials{cond_ix}(ch_ix,:,:)),...
%                 'Color',cond_colors{cond_ix},'LineWidth',plt.butterfly_width,...
%                 'LineStyle',cond_styles{cond_ix});
%         end
%     end
% 
    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(rois{1}.time{1}, plt_means(cond_ix,:), sems(cond_ix,:),...
            {'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    ylims = ylim;
%     if strcmp(plt.sig_type,'line')
%         data_lim = [min(min(means-sems)) max(max(means+sems))];
%     end
    
    % Plot Extra Features (events, significance)
    for evnt_ix = 1:numel(an.event_type)
        main_lines(numel(cond_lab)+evnt_ix) = line([0 0],ylim,...
            'LineWidth',plt.evnt_width(evnt_ix),'Color',plt.evnt_color{evnt_ix},'LineStyle',plt.evnt_style{evnt_ix});
    end
    
%     if strcmp(plt.sig_type,'line')
%         leg_lab = [cond_lab an.event_type grp_lab(sig_grp)];
%     else
        leg_lab = [cond_lab an.event_type];
%     end
    
    % Axes and Labels
%     ax.YLim          = ylims; %!!! change for plt.sigType=line
    ax.YLabel.String = 'uV';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    ax.Title.String  = w2s{s}.label{ch_ix};
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
%     fig_name = ['GRP_' conditions '_' an_id '_' roi_erps{1,1}.label{ch_ix}];
% %     [plot_rc,~] = fn_num_subplots(numel(w2.label));
% %     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     
%     figure('Name',fig_name,'units','normalized',...
%         'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
%     plot_info.fig        = gcf;
%     plot_info.x_step     = plt.x_step_sz*proc.resample_freq;
%     plot_info.x_lab      = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
%     plot_info.y_lab      = 'uV';
%     plot_info.legend_loc = plt.legend_loc;
%     plot_info.sig_alpha  = plt.sig_alpha;
%     plot_info.sig_color  = plt.sig_color;
%     % Stimulus plotting params
%     event_info.time      = -plt.plt_lim(1)*proc.resample_freq;
%     event_info.name      = {an.event_type};
%     event_info.width     = plt.evnt_width;
%     event_info.color     = {plt.evnt_color};
%     event_info.style     = {plt.evnt_style};
%     % Condition plotting params
%     cond_info.name       = cond_lab;
%     cond_info.style      = cond_style;
%     cond_info.color      = cond_colors;
%     cond_info.alpha      = repmat(plt.errbar_alpha,[1 numel(cond_lab)]);
%     
% %     subplot(plot_rc(1),plot_rc(2),ch_ix);
%     plot_info.ax     = gca;
%     plot_info.title  = strcat(roi_erps{1,1}.label{ch_ix});
%     plot_info.legend = plt.legend;
%     
%     % Compute means and variance
%     means = NaN([numel(cond_lab) size(roi_erps{1,1}.avg,2)]);
%     sem = NaN([numel(cond_lab) size(roi_erps{1,1}.avg,2)]);
%     for cond_ix = 1:numel(cond_lab)
%         sbj_means = NaN([numel(SBJs) size(means,2)]);
%         for s = 1:numel(SBJs)
%             sbj_means(s,:) = roi_erps{s,cond_ix}.avg(ch_ix,:);
%         end
%         means(cond_ix,:) = mean(sbj_means,1);
%         sem(cond_ix,:) = std(sbj_means)./sqrt(numel(SBJs));
%     end
%     % Find significant time periods
% %     if sum(w2.mask(ch_ix,:))>0
% %         sig_ch(ch_ix) = 1;
% %         mask_chunks = fn_find_chunks(w2.mask(ch_ix,:));
% %         sig_chunks = mask_chunks;
% %         sig_chunks(w2.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
% %         % If w2 and roi_erps aren't on same time axis, adjust sig_chunk indices
% %         if (size(w2.time,2)~=size(roi_erps{1}.time,2)) || (sum(stat.time==roi_erps{1}.time)~=numel(stat.time))
% %             for chunk_ix = 1:size(sig_chunks,1)
% %                 sig_chunks(chunk_ix,1) = find(roi_erps{1}.time==stat.time(sig_chunks(chunk_ix,1)));
% %                 sig_chunks(chunk_ix,2) = find(roi_erps{1}.time==stat.time(sig_chunks(chunk_ix,2)));
% %             end
% %         end
% %         fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
% %                                                                 stat.label{ch_ix},size(sig_chunks,1));
% %         fn_plot_ts_errbr_sig(plot_info,means,sem,sig_chunks,event_info,cond_info);
% %     else
% %         fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat.label{ch_ix});
%     fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
% %     end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_filename]);
    end
end

%% Save out list of channels with significant differences
% if write_report
%     sig_report_fname = [fig_dir 'ch_sig_list.csv'];
%     sig_report = fopen(sig_report_fname,'w');
%     fprintf(sig_report,'%s\n',an_id);
%     fprintf(sig_report,'label,%s\n',conditions);
%     for ch_ix = 1:numel(stat.label)
%         fprintf(sig_report,'%s,%.0f\n',stat.label{ch_ix},sig_ch(ch_ix));
%     end
%     fclose(sig_report);
% end

end
