function SBJ09b_TFR_plot_stats(SBJ,conditions,proc_id,an_id,plt_id,save_fig,fig_vis)
% Plots TFRs computed in SBJ09a_TFR_stats
%   3 TFRs in a row: (1) con; (2) inc; (3) inc-con with stats
% clear all; %close all;
error('not adapted to stat_id yet');

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

%% Data Preparation
% Load Data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

stats_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id,'.mat');
load(stats_filename);

% Select Conditions of Interest
[cond_lab, ~, ~] = fn_condition_label_styles(conditions);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
% Compute mean RT per condition
for cond_ix = 1:numel(cond_lab)
    if strcmp(event_type,'stim')
        RT_means{cond_ix} = mean(trial_info.response_time(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1));
    else
        RT_means{cond_ix} = 0;
    end
end

% Load ROI and GM/WM info
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_Dx.mat'];
load(elec_fname);
% Select TFR elecs
cfgs = []; cfgs.channel = stat.label;
elec = fn_select_elec(cfgs, elec);
% Get ROI for figure label
elec.roi = fn_atlas2roi_labels(elec.atlas_label,'Dx','ROI');

%% Prep Data for Plotting
% Obtain the descriptive stats
tfr_avg = {};
cfg_avg = [];
cfg_avg.channel = 'all';
cfg_avg.latency = trial_lim_s;
for cond_ix = 1:numel(cond_lab)
    % Trim to plotting epoch
    tfr_avg{cond_ix} = ft_freqdescriptives(cfg_avg,tfr{cond_ix});
end

% Take the difference
tfr_diff = tfr_avg{1};
% for cond_ix = 2:numel(cond_lab)
tfr_diff.powspctrm = tfr_avg{2}.powspctrm - tfr_avg{1}.powspctrm;
% end
cfg_trim = [];
cfg_trim.latency = stat_lim;
tfr_diff = ft_selectdata(cfg_trim,tfr_diff);
tfr_diff.statmask = stat.mask;
% stat.raweffect = tfr_diff{1}; % I think this was for ft_clusterplot, maybe only for topos?

%% Plot results
fig_dir = [root_dir 'PRJ_Stroop/results/TFR/' SBJ '/' conditions '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure per channel
sig_ch = {};
for ch_ix = 1:numel(stat.label)
    fig_name = [SBJ '_' conditions '_' an_id '_' ...
        elec.roi{strcmp(elec.label,tfr{cond_ix}.label{ch_ix})} '_' tfr{cond_ix}.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.5],'Visible',fig_vis);
    
    cfg_s = [];
    cfg_s.channel = tfr{cond_ix}.label{ch_ix};
%     cfg_lay = [];
%     cfg_lay.layout = 'vertical';

    % Condition Plots
    for cond_ix = 1:length(cond_lab)
        subplot(1,numel(cond_lab)+1,cond_ix);   %!!! assumes only 2 conditions
        tfr_plot = ft_selectdata(cfg_s, tfr_avg{cond_ix});
%         layout = ft_prepare_layout(cfg_lay,tfr_plot);
        freqs = tfr_plot.freq;
        tfr_plot.freq = 1:1:numel(freqs);
        
        cfgraw.latency  = trial_lim_s;
        cfgraw.baseline = bsln_lim; %should be in sec
        ft_singleplotTFR(cfgraw, tfr_plot);
        line([RT_means{cond_ix} RT_means{cond_ix}],ylim,'Color','k','LineStyle','--','LineWidth',2);
        
        ax = gca;
        ax.YTick = 1:tick_step:numel(freqs);
        ax.YTickLabel = round(freqs(1:tick_step:end));
        title(cond_lab{cond_ix});
    end
    
    % Difference Plot with Stats
    subplot(1,numel(cond_lab)+1,numel(cond_lab)+1);
    freqs = tfr_diff.freq;
    tfr_diff.freq        = 1:1:numel(freqs);
    cfgdif.latency       = stat_lim;
    cfgdif.channel       = tfr_diff.label{ch_ix};
    ft_singleplotTFR(cfgdif,tfr_diff);
    
    ax = gca;
    ax.YTick = 1:tick_step:numel(freqs);
    ax.YTickLabel = round(freqs(1:tick_step:end));
    title([cond_lab{2} '-' cond_lab{1}]);
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
    if sum(squeeze(stat.mask(ch_ix,:,:)))>0
        sig_ch = {sig_ch{:} stat.label{ch_ix}};
        % Link this figure to sig dir
        sig_dir = [fig_dir 'sig_ch/'];
        if ~exist(sig_dir,'dir')
            mkdir(sig_dir);
        end
        cd(sig_dir);
        link_cmd = ['ln -s ../' fig_name ' .'];
        system(link_cmd);
    end
end

% Save out list of channels with significant differences
sig_report_filename = [fig_dir 'sig_ch_list.txt'];
sig_report = fopen(sig_report_filename,'a');
fprintf(sig_report,'%s\n',an_id,an_id);
fprintf(sig_report,'%s\n',sig_ch{:});
fclose(sig_report);
fprintf('===================================================\n');
fprintf('--- SIGNIFICANT CHANNELS: %s ------------------\n',sig_ch{:});
fprintf('===================================================\n');

end

%% %% Cluster plots
% cfg_clust = [];
% cfg_clust.alpha  = 0.025;
% cfg_clust.parameter = 'raweffect';
% cfg_clust.zlim   = [-1e-27 1e-27];
% cfg_clust.layout = layout;
% ft_clusterplot(cfg_clust, stat);

%% cfg_s = [];
% cfg_s.latency = [-0.5 2];
% tfr_plot = ft_selectdata(cfg_s, tfr);
% figure;
% for ch_ix = 1:length(tfr.label)
%     subplot(2,1,ch_ix);
%     cfg = [];
%     cfg.baseline = [-.2 0];%bsln_lim; %should be in sec
%     cfg.baselinetype = 'db';
%     cfg.channel = tfr.label{ch_ix};
%     cfg.colorbar = 'yes';
%     cfg.zlim = 'maxabs'; %[-5 5];
%     cfg.parameter = 'powspctrm';
%     % cfg.maskparameter = 'trial';
%     % cfg.maskstyle = 'outline';
%     cfg.showlabels   = 'yes';
%     ft_singleplotTFR(cfg, tfr_plot);
% end
