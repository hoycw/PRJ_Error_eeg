function SBJ07b_ERP_plot_SR_stats_stack(SBJ,conditions,an_id_s,an_id_r,plt_id,save_fig,fig_vis)
% Plots ERPs computed in SBJ07a_ERP_stats
% clear all; %close all;
plt_vars.clim_perc = [5 95];

%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
% eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

event_lab = {'stim', 'resp'};

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n));
end

[cond_mat,good_trials] = find(cond_idx);
cond_mat = horzcat(cond_mat,round(1000*trial_info.response_time(good_trials)),[1:numel(cond_mat)]');
cond_mat = sortrows(cond_mat,[1 2]);

% Load Data
stats_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_s,'.mat');
stats_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_r,'.mat');
tmp = load(stats_filename1,'stat'); stat{1} = tmp.stat;
tmp = load(stats_filename2,'stat'); stat{2} = tmp.stat;
tmp = load(stats_filename1,'roi_erp'); erp{1,1} = tmp.roi_erp{1}; erp{1,2} = tmp.roi_erp{2};
tmp = load(stats_filename2,'roi_erp'); erp{2,1} = tmp.roi_erp{1}; erp{2,2} = tmp.roi_erp{2};
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(erp{1,1}.time)-1)/(erp{1,1}.time(end)-erp{1,1}.time(1));
if ~isempty(setdiff(stat{1}.label,stat{2}.label))
    error('ERROR: channels do not match between the two analyses!');
end

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
erp{1,1} = ft_selectdata(cfg_trim,erp{1,1});
erp{1,2} = ft_selectdata(cfg_trim,erp{1,2});
cfg_trim.latency = plt_vars.plt_lim_R;
erp{2,1} = ft_selectdata(cfg_trim,erp{2,1});
erp{2,2} = ft_selectdata(cfg_trim,erp{2,2});

% Compute mean RT per condition
RTs = round(1000*trial_info.response_time); % converts sec to ms
for cond_ix = 1:numel(cond_lab)
    RT_means{cond_ix} = mean(RTs(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1));
    % Add in the baseline offset to plot correctly
    RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*1000;
end

%% Plot Results
% NO! I will plot whatever channels I ran the stats on (why else did I run them?)
% % Select data to plot only ROI channels
% cfgs = [];
% cfgs.channel = SBJ_vars.ch_lab.ROI;
% stat = ft_selectdata(cfgs,stat);
% for an_ix = 1:numel(cond_lab)
%     erp{an_ix} = ft_selectdata(cfgs,erp{an_ix});
% end

fig_dir = [root_dir 'PRJ_Stroop/results/ERP/' SBJ '/' conditions '/SR/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = {};
for ch_ix = 1:numel(stat{1}.label)
    % Plot parameters
%     probe_name = stat.label{ch_ix}(regexp(stat.label{ch_ix},'\D'));
%     probe_name(strfind(probe_name,'-')) = [];
    
    fig_name = [SBJ '_' conditions '_SR_' stat{1}.label{ch_ix}];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.5],'Visible',fig_vis);   %this size is for single plots
    plot_info.fig        = gcf;
    plot_info.x_step     = plt_vars.x_step_sz*sample_rate;
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);
    
    % Plot Mean Traces
    for sr_ix = 1:2
        subplot(2,2,sr_ix);
        plot_info.ax     = gca;
        plot_info.title  = [stat{sr_ix}.label{ch_ix} ':' event_lab{sr_ix}];
        plot_info.legend = plt_vars.legend;
        if strcmp(event_lab{sr_ix},'stim')
            plot_info.x_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{sr_ix}, cond_lab{:}};
            event_info.color = {[0 0 0], cond_colors{:}};
            event_info.width = repmat(plt_vars.evnt_width,[1 numel(event_info.name)]);
            event_info.style = repmat({plt_vars.evnt_style},[1 numel(event_info.name)]);
            event_info.time  = [-plt_vars.plt_lim_S(1)*sample_rate, RT_means{:}];
        else
            plot_info.x_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Stimulus plotting params
            event_info.name  = {event_lab{sr_ix}};
            event_info.width = plt_vars.evnt_width;
            event_info.color = {plt_vars.evnt_color};
            event_info.style = {plt_vars.evnt_style};
            event_info.time  = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        
        % Compute means and variance
        means = NaN([numel(cond_lab) size(erp{sr_ix,1}.time,2)]);
        var = NaN([numel(cond_lab) size(erp{sr_ix,1}.time,2)]);
        for an_ix = 1:numel(cond_lab)
            means(an_ix,:) = squeeze(mean(erp{sr_ix,an_ix}.trial(:,ch_ix,:),1));
            var(an_ix,:) = squeeze(std(erp{sr_ix,an_ix}.trial(:,ch_ix,:),[],1)./...
                                                sqrt(size(erp{sr_ix,an_ix}.trial,1)))';
        end
        % Find significant time periods
        if sum(stat{sr_ix}.mask(ch_ix,:))>0
            sig_ch = {sig_ch{:} stat{sr_ix}.label{ch_ix}};
            mask_chunks = fn_find_chunks(stat{sr_ix}.mask(ch_ix,:));
            sig_chunks = mask_chunks;
            sig_chunks(stat{sr_ix}.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
            % If stat and erp aren't on same time axis, adjust sig_chunk indices
            if (size(stat{sr_ix}.time,2)~=size(erp{sr_ix,1}.time,2)) || ...
                                (sum(stat{sr_ix}.time==erp{sr_ix,1}.time)~=numel(stat{sr_ix}.time))
                for chunk_ix = 1:size(sig_chunks,1)
                    sig_chunks(chunk_ix,1) = find(erp{sr_ix,1}.time==stat{sr_ix}.time(sig_chunks(chunk_ix,1)));
                    sig_chunks(chunk_ix,2) = find(erp{sr_ix,1}.time==stat{sr_ix}.time(sig_chunks(chunk_ix,2)));
                end
            end
            fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                stat{sr_ix}.label{ch_ix},size(sig_chunks,1));
            fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
        else
            fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat{sr_ix}.label{ch_ix});
            fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
        end
    end
    
    % Plot Single Trial Stacks
    % Get color limits
    clims = NaN([2 2]);
    for sr_ix = 1:2
        % mins across conditions
        clims(sr_ix,1) = prctile(reshape(erp{sr_ix,1}.trial(:,ch_ix,:),...
            [1 sum(cond_mat(:,1)==1)*numel(erp{sr_ix,1}.time)]),plt_vars.clim_perc(1));
        for cond_ix = 2:numel(cond_lab)
            clims(sr_ix,1) = min(clims(sr_ix,1),prctile(reshape(erp{sr_ix,2}.trial(:,ch_ix,:),...
                [1 sum(cond_mat(:,1)==cond_ix)*numel(erp{sr_ix,2}.time)]),plt_vars.clim_perc(1)));
        end
        % maxes across conditions
        clims(sr_ix,2) = prctile(reshape(erp{sr_ix,1}.trial(:,ch_ix,:),...
            [1 sum(cond_mat(:,1)==1)*numel(erp{sr_ix,1}.time)]),plt_vars.clim_perc(2));
        for cond_ix = 2:numel(cond_lab)
            clims(sr_ix,2) = max(clims(sr_ix,2),prctile(reshape(erp{sr_ix,2}.trial(:,ch_ix,:),...
                [1 sum(cond_mat(:,1)==cond_ix)*numel(erp{sr_ix,2}.time)]),plt_vars.clim_perc(2)));
        end
    end
    clims = [min(clims(:,1)) max(clims(:,2))];
    for sr_ix = 1:2
        subplot(2,2,2+sr_ix);
        hold on;
        
        % Plot Single Trials Per Condition
        imagesc(vertcat(squeeze(erp{sr_ix,1}.trial(:,ch_ix,:)),squeeze(erp{sr_ix,2}.trial(:,ch_ix,:))));
        set(gca,'YDir','normal');
        if strcmp(event_lab{sr_ix},'stim')
            x_tick_lab = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            event_time = -plt_vars.plt_lim_S(1)*sample_rate;
            y_ix = 1;
            for cond_ix = 1:numel(cond_lab)
                idx = cond_mat(:,1)==cond_ix;
                scat = scatter(cond_mat(idx,2)-plt_vars.plt_lim_S(1)*sample_rate,find(idx),...
                    'MarkerFaceColor',[cond_colors{cond_ix}],'MarkerEdgeColor','k');
            end
        else
            x_tick_lab = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            event_time = -plt_vars.plt_lim_R(1)*sample_rate;
        end
        ylim([1 size(cond_mat,1)]);
        event_line = line([event_time event_time],ylim,...
            'LineWidth',plt_vars.evnt_width,'Color','k');

        % Plotting parameters
        ax = gca;
%         ax.legend = plt_vars.legend;
        ax.Title.String  = [erp{sr_ix,1}.label{ch_ix} ' trials'];
        ax.XLim          = [0,numel(erp{sr_ix,1}.time)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:numel(erp{sr_ix,1}.time);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Trials';
%         legend([roi_lines{logical(roi_flag)}],lgd_lab{logical(roi_flag)},'Location',lgd_loc);
        cbar = colorbar;
        caxis(clims);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

%% Save out list of channels with significant differences
sig_report_filename = [fig_dir 'sig_ch_list.txt'];
sig_report = fopen(sig_report_filename,'a');
fprintf(sig_report,'%s - %s\n',an_id_s,an_id_r);
fprintf(sig_report,'%s\n',sig_ch{:});
fclose(sig_report);

end