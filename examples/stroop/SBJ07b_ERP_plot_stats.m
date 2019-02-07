function SBJ07b_ERP_plot_stats(SBJ,conditions,an_id_s,an_id_r,plt_id,save_fig,fig_vis)
% Plots ERPs computed in SBJ07a_ERP_stats
% clear all; %close all;

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);

stats_filename = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_ROI_',an_id,'.mat');
load(stats_filename);

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(roi_erp{1}.time)-1)/(roi_erp{1}.time(end)-roi_erp{1}.time(1));

% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_lim;
roi_erp{1} = ft_selectdata(cfg_trim,roi_erp{1});
roi_erp{2} = ft_selectdata(cfg_trim,roi_erp{2});

%% Plot Results
% NO! I will plot whatever channels I ran the stats on (why else did I run them?)
% % Select data to plot only ROI channels
% cfgs = [];
% cfgs.channel = SBJ_vars.ch_lab.ROI;
% stat = ft_selectdata(cfgs,stat);
% for an_ix = 1:numel(cond_lab)
%     roi_erp{an_ix} = ft_selectdata(cfgs,roi_erp{an_ix});
% end

fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/ERP/' SBJ '/' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = {};
for ch_ix = 1:numel(stat.label)
    % Plot parameters
%     probe_name = stat.label{ch_ix}(regexp(stat.label{ch_ix},'\D'));
%     probe_name(strfind(probe_name,'-')) = [];
    
    fig_name = [SBJ '_' conditions '_' an_id '_' stat.label{ch_ix}];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
    plot_info.fig        = gcf;
    plot_info.x_step     = plt_vars.x_step_sz*sample_rate;
    plot_info.x_lab      = plt_lim(1):plt_vars.x_step_sz:plt_lim(2);
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    % Stimulus plotting params
    event_info.time      = -plt_lim(1)*sample_rate;
    event_info.name      = event_type;
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = plt_vars.evnt_color;
    event_info.style     = plt_vars.evnt_style;
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);
    
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    plot_info.ax     = gca;
    plot_info.title  = stat.label{ch_ix};
    plot_info.legend = plt_vars.legend;
    
    % Compute means and variance
    means = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
    var = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
    for an_ix = 1:numel(cond_lab)
        means(an_ix,:) = roi_erp{an_ix}.avg(ch_ix,:);
        var(an_ix,:) = squeeze(std(roi_erp{an_ix}.trial(:,ch_ix,:),[],1)./sqrt(size(roi_erp{an_ix}.trial,1)))';
    end
    % Find significant time periods
    if sum(stat.mask(ch_ix,:))>0
        sig_ch = {sig_ch{:} stat.label{ch_ix}};
        mask_chunks = fn_find_chunks(stat.mask(ch_ix,:));
        sig_chunks = mask_chunks;
        sig_chunks(stat.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
        % If stat and roi_erp aren't on same time axis, adjust sig_chunk indices
        if (size(stat.time,2)~=size(roi_erp{1}.time,2)) || (sum(stat.time==roi_erp{1}.time)~=numel(stat.time))
            for chunk_ix = 1:size(sig_chunks,1)
                sig_chunks(chunk_ix,1) = find(roi_erp{1}.time==stat.time(sig_chunks(chunk_ix,1)));
                sig_chunks(chunk_ix,2) = find(roi_erp{1}.time==stat.time(sig_chunks(chunk_ix,2)));
            end
        end
        fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                                                                stat.label{ch_ix},size(sig_chunks,1));
        fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
    else
        fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat.label{ch_ix});
        fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

% Save out list of channels with significant differences
sig_report_filename = [fig_dir 'sig_ch_list.txt'];
sig_report = fopen(sig_report_filename,'a');
fprintf(sig_report,'%s\n',an_id);
fprintf(sig_report,'%s\n',sig_ch{:});
fclose(sig_report);

end
