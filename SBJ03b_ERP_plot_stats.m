function SBJ03b_ERP_plot_stats(SBJ,conditions,pipeline_id,an_id,plt_id,save_fig,fig_vis)
% clear all; %close all;

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';app_dir='/Users/SCS22/Documents/MATLAB/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

[cond_lab, cond_colors, cond_style, ~] = fn_condition_label_styles(conditions);

stats_filename = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_',an_id,'.mat');
load(stats_filename);

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(roi_erp{1}.time)-1)/(roi_erp{1}.time(end)-roi_erp{1}.time(1));

cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim;
for an_ix = 1:numel(roi_erp)
    % Get rid of .trial field to avoid following bug that tosses .avg, .var, and .dof:
    %   "Warning: timelock structure contains field with and without repetitions"
    roi_erp{an_ix} = rmfield(roi_erp{an_ix},'trial');    
    
    % Trim data to plotting epoch
    roi_erp{an_ix} = ft_selectdata(cfg_trim,roi_erp{an_ix});    
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' SBJ '/' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = zeros(size(stat.label));
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
    plot_info.x_lab      = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    plot_info.legend_loc = plt_vars.legend_loc;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.sig_color  = plt_vars.sig_color;
    % Stimulus plotting params
    event_info.time      = -plt_vars.plt_lim(1)*sample_rate;
    event_info.name      = {event_type};
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = {plt_vars.evnt_color};
    event_info.style     = {plt_vars.evnt_style};
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);
    
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    plot_info.ax     = gca;
    plot_info.title  = strcat(stat.label{ch_ix});
    plot_info.legend = plt_vars.legend;
    
    % Compute means and variance
    means = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
    sem = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
    for an_ix = 1:numel(cond_lab)
        means(an_ix,:) = roi_erp{an_ix}.avg(ch_ix,:);
        sem(an_ix,:) = squeeze(sqrt(roi_erp{an_ix}.var(ch_ix,:))./sqrt(numel(roi_erp{an_ix}.cfg.previous.trials)))';
    end
    % Find significant time periods
    if sum(stat.mask(ch_ix,:))>0
        sig_ch(ch_ix) = 1;
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
        fn_plot_ts_errbr_sig(plot_info,means,sem,sig_chunks,event_info,cond_info);
    else
        fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat.label{ch_ix});
        fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
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
sig_report_filename = [fig_dir 'ch_sig_list.csv'];
sig_report = fopen(sig_report_filename,'w');
fprintf(sig_report,'%s\n',an_id);
fprintf(sig_report,'label,%s\n',conditions);
for ch_ix = 1:numel(stat.label)
    fprintf(sig_report,'%s,%.0f\n',stat.label{ch_ix},sig_ch(ch_ix));
end
fclose(sig_report);

end