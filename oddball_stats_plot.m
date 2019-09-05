function oddball_stats_plot(SBJ, proc_id, plt_id, an_id, fig_vis)
%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/', ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip_private']);
addpath(ft_dir);
ft_defaults

%% Load preprocessed data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_oddball_' proc_id '_clean.mat'];
load(clean_bhv_fname);
data_cleanname = [SBJ_vars.dirs.preproc SBJ '_clean_oddball_' proc_id '.mat'];
load(data_cleanname);

%find index of conditions of interest
std_cdx = find(clean_oddball.trialinfo == 1);
tar_cdx = find(clean_oddball.trialinfo == 2);
odd_cdx = find(clean_oddball.trialinfo == 3);

% Select Channel(s)
cfgs = [];
cfgs.channel = ROI;
roi = ft_selectdata(cfgs, clean_oddball);

roi_erp = {};
cfgavg = [];
cfgavg.keeptrials = 'yes';
cfgavg.trials = std_cdx;
roi_erp{1} = ft_timelockanalysis(cfgavg,roi);
ntrials(1) = size(roi_erp{1}.trial, 1);
cfgavg.trials = tar_cdx;
roi_erp{2} = ft_timelockanalysis(cfgavg,roi);
ntrials(2) = size(roi_erp{2}.trial, 1);
cfgavg.trials = odd_cdx;
roi_erp{3} = ft_timelockanalysis(cfgavg,roi);
ntrials(3) = size(roi_erp{3}.trial, 1);
%make design matrix and stats

design = zeros(2,ntrials(2) + ntrials(3));
for cond_ix = 2:3
    if cond_ix==2
        design(1,1:ntrials(cond_ix)) = cond_ix;                                % Conditions (Independent Variable)
        design(2,1:ntrials(cond_ix)) = 1:ntrials(cond_ix);                    % Trial Numbers
    else
        design(1,(ntrials(cond_ix-1)+1):sum(ntrials(cond_ix-1:cond_ix)))= cond_ix; % Conditions (Independent Variable)
        design(2,(ntrials(cond_ix-1)+1):sum(ntrials(cond_ix-1:cond_ix)))= 1:ntrials(cond_ix);
    end
end

% Calculate statistics

cfg_stat.design           = design;
[stat] = ft_timelockstatistics(cfg_stat, roi_erp{2}, roi_erp{3});
%{
%% Stats
tar = squeeze(mean(roi_erp{2}.trial));
odd = squeeze(mean(roi_erp{3}.trial));
%One Sided
[h, p, ci, stats] = ttest(tar-odd, 0, 0.05)
%Two Sided
[h,p,ci,stats] = ttest2(tar, odd)
%%Plot
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/oddball' SBJ  '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
%}
for ch_ix = 1:numel(ROI)
    fig_name = [SBJ '_oddball_' ROI{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible', 'on');   %this size is for single plots
    plot_info.fig = gcf;
    
    plot(squeeze(mean(roi_erp{1}.trial)));
    hold on
    plot(squeeze(mean(roi_erp{2}.trial)));
    hold on
    plot(squeeze(mean(roi_erp{3}.trial)));
    plot_info.ax     = gca;
    legend('standard', 'target', 'odd');
end


for ch_ix = 1:numel(ROI)  %numel roi if channel
% Stimulus plotting params
%event_info.time      = -plt_vars.plt_lim(1)*sample_rate;
%{
%%Plot info that came from original fucntion -- req'd for error bar
funciton
event_info.name      = {event_type};
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = {plt_vars.evnt_color};
    event_info.style     = {plt_vars.evnt_style};
    % Condition plotting params
    cond_info.name       = ['standard', 'target', 'odd'];
    cond_info.style      = {'-', '-', '-'};
    cond_info.color      = {'r', 'g', 'b'};
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 3]);
    
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    plot_info.ax     = gca;
    %plot_info.title  = strcat(stat.label{ch_ix});
    plot_info.legend = plt_vars.legend;
    %}
    % Compute means and variance
    for cond_ix = 1:3
         roi_erp{cond_ix}.avg = squeeze(mean(roi_erp{cond_ix}.trial));
         roi_erp{cond_ix}.var = squeeze(var(roi_erp{cond_ix}.trial));
    end
    means = NaN([3 size(roi_erp{1}.avg,2)]);
    sem = NaN([3 size(roi_erp{1}.avg,2)]);
    for an_ix = 1:3
        means(an_ix,:) = roi_erp{an_ix}.avg(ch_ix,:);
        sem(an_ix,:) = squeeze(sqrt(roi_erp{an_ix}.var(ch_ix,:))./sqrt(numel(roi_erp{an_ix}.cfg.previous.trials)))';
    end  
 
  %Find significant time periods
    if sum(stat.mask(ch_ix,:))>0
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
        %fn_plot_ts_errbr_sig(plot_info,means,sem,sig_chunks,event_info,cond_info);
    else
        fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat.label{ch_ix});
        %fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
    end
%% Save stats variable
stats_fname = [SBJ_vars.dirs.proc SBJ '_stats_' an_id '.mat'];
save(stats_fname, '-v7.3', 'stat');

% Save figure    
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/oddball' SBJ  '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
fig_filename = [fig_dir fig_name '.png'];
fprintf('Saving %s\n',fig_filename);
saveas(gcf,fig_filename);
%eval(['export_fig ' fig_filename]);

end

