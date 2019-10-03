function SBJ03c_ERP_plot_group(SBJs,conditions,proc_id,an_id,plt_id,save_fig,varargin)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';app_dir='/Users/SCS22/Documents/MATLAB/';
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

[cond_lab, cond_colors, cond_style, ~] = fn_condition_label_styles(conditions);

% Load Data
cfg_trim = [];
cfg_trim.latency = plt.plt_lim;
SBJs_vars = cell(size(SBJs));
roi_erps  = cell([numel(SBJs) numel(cond_lab)]);
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    SBJs_vars{s} = SBJ_vars;
    
    tmp = load([SBJs_vars{s}.dirs.SBJ,'04_proc/',SBJs{s},'_',conditions,'_',an_id,'.mat']);
    for cond_ix = 1:numel(cond_lab)
        roi_erps{s,cond_ix} = ft_selectdata(cfg_trim,tmp.roi_erp{cond_ix});
    end
    clear SBJ_vars
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/GRP/' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
sig_ch = zeros(size(stat.label));
for ch_ix = 1:numel(stat.label)
    fig_name = ['GRP_' conditions '_' an_id '_' roi_erps{1,1}.label{ch_ix}];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
    plot_info.fig        = gcf;
    plot_info.x_step     = plt.x_step_sz*proc.resample_freq;
    plot_info.x_lab      = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    plot_info.y_lab      = 'uV';
    plot_info.legend_loc = plt.legend_loc;
    plot_info.sig_alpha  = plt.sig_alpha;
    plot_info.sig_color  = plt.sig_color;
    % Stimulus plotting params
    event_info.time      = -plt.plt_lim(1)*proc.resample_freq;
    event_info.name      = {an.event_type};
    event_info.width     = plt.evnt_width;
    event_info.color     = {plt.evnt_color};
    event_info.style     = {plt.evnt_style};
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = cond_style;
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt.errbar_alpha,[1 numel(cond_lab)]);
    
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    plot_info.ax     = gca;
    plot_info.title  = strcat(roi_erps{1,1}.label{ch_ix});
    plot_info.legend = plt.legend;
    
    % Compute means and variance
    means = NaN([numel(cond_lab) size(roi_erps{1,1}.avg,2)]);
    sem = NaN([numel(cond_lab) size(roi_erps{1,1}.avg,2)]);
    for cond_ix = 1:numel(cond_lab)
        sbj_means = NaN([numel(SBJs) size(means,2)]);
        for s = 1:numel(SBJs)
            sbj_means(s,:) = roi_erps{s,cond_ix}.avg(ch_ix,:);
        end
        means(cond_ix,:) = mean(sbj_means,1);
        sem(cond_ix,:) = std(sbj_means)./sqrt(numel(SBJs));
    end
    % Find significant time periods
%     if sum(stat.mask(ch_ix,:))>0
%         sig_ch(ch_ix) = 1;
%         mask_chunks = fn_find_chunks(stat.mask(ch_ix,:));
%         sig_chunks = mask_chunks;
%         sig_chunks(stat.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
%         % If stat and roi_erps aren't on same time axis, adjust sig_chunk indices
%         if (size(stat.time,2)~=size(roi_erps{1}.time,2)) || (sum(stat.time==roi_erps{1}.time)~=numel(stat.time))
%             for chunk_ix = 1:size(sig_chunks,1)
%                 sig_chunks(chunk_ix,1) = find(roi_erps{1}.time==stat.time(sig_chunks(chunk_ix,1)));
%                 sig_chunks(chunk_ix,2) = find(roi_erps{1}.time==stat.time(sig_chunks(chunk_ix,2)));
%             end
%         end
%         fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
%                                                                 stat.label{ch_ix},size(sig_chunks,1));
%         fn_plot_ts_errbr_sig(plot_info,means,sem,sig_chunks,event_info,cond_info);
%     else
%         fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat.label{ch_ix});
    fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
%     end
    
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
