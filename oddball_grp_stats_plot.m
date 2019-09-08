function oddball_grp_stats_plot(SBJ, proc_id, plt_id, an_id, fig_vis, save_fig, fig_ftype)

% load all SBJ roi_erp
roi_erp = cell(size(SBJs));
for s = 1:numel(SBJs)
    tmp = load(stat_fname); roi_erp{s} = tmp.roi_erp;
end
cfg_avg = [];
cfg_avg.channel = {'all'};
cfg_avg.latency = an.stat_lim;
cfg_avg.keepindividual = 'yes';
cfg_avg.method = 'across';
cfg_avg.parameter = 'avg';
gstat = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    gstat = ft_timelockgrandaverage(cfg_avg, roi_erp{:}{cond_ix});
end

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);
addpath(ft_dir);
ft_defaults

%% Load preprocessed data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat'];
load(clean_bhv_fname);
data_cleanname = [SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat'];
load(data_cleanname);

% Condition info
cond_lab = {'std','tar','odd'};
cond_colors = {[0 0 1],[0 1 0],[1 0 0]};
cond_comp = {[1 2], [1 3]}; % std vs. tar; std vs. odd

%% Compute ERPs
% Select Channel(s)
cfgs = [];
cfgs.channel = an.ROI;
roi = ft_selectdata(cfgs, clean_trials);

roi_erp  = cell(size(cond_lab));
n_trials = zeros(size(cond_lab));
cfgavg = [];
cfgavg.keeptrials = 'yes';
for cond_ix = 1:numel(cond_lab)
    cfgavg.trials = find(clean_trials.trialinfo == cond_ix);
    roi_erp{cond_ix} = ft_timelockanalysis(cfgavg,roi);
    roi_erp{cond_ix}.avg = squeeze(mean(roi_erp{cond_ix}.trial,1))';
    roi_erp{cond_ix}.var = squeeze(var(roi_erp{cond_ix}.trial,[],1))';
    n_trials(cond_ix) = size(roi_erp{cond_ix}.trial, 1);
end

%% Calculate statistics
stat   = cell(size(cond_comp));
design = cell(size(cond_comp));
for st_ix = 1:numel(cond_comp)
    cond_ixs = cond_comp{st_ix};
    
    % Create make design matrix and stats
    design{st_ix} = zeros(2, n_trials(cond_ixs(1)) + n_trials(cond_ixs(2)));
    for c_ix = 1:2
        if c_ix==1
            design{st_ix}(1,1:n_trials(cond_ixs(c_ix))) = cond_ixs(c_ix);                               % Conditions (Independent Variable)
            design{st_ix}(2,1:n_trials(cond_ixs(c_ix))) = 1:n_trials(cond_ixs(c_ix));                   % Trial Numbers
        else
            design{st_ix}(1,n_trials(cond_ixs(c_ix-1))+1:sum(n_trials(cond_ixs))) = cond_ixs(c_ix); % Conditions (Independent Variable)
            design{st_ix}(2,n_trials(cond_ixs(c_ix-1))+1:sum(n_trials(cond_ixs)))= 1:n_trials(cond_ixs(c_ix));
        end
    end
    
    % Compute stats between ERPs
    cfg_stat.design           = design{st_ix};
    [stat{st_ix}] = ft_timelockstatistics(cfg_stat, roi_erp{cond_ixs(1)}, roi_erp{cond_ixs(2)});
end
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

%% Plot ERPs and stats
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' SBJ '/oddball/' an_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
end

for ch_ix = 1:numel(an.ROI)
    fig_name = [SBJ '_oddball_' an.ROI{ch_ix}];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible', fig_vis);   %this size is for single plots
    
    % General plotting params
    plot_info.fig        = f;
    plot_info.sig_color  = plt_vars.sig_color;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.x_step = plt_vars.x_step_sz*roi.fsample;
    plot_info.x_lab  = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    plot_info.y_lab  = 'uV';
    plot_info.legend = plt_vars.legend;
    plot_info.legend_loc = plt_vars.legend_loc;
    % Event plotting params
    event_info.name      = {an.event_type};
    [~,event_info.time]  = min(abs(roi_erp{1}.time-0));
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = {plt_vars.evnt_color};
    event_info.style     = {plt_vars.evnt_style};
    % Condition plotting params
    cond_info.name       = cond_lab;
    cond_info.style      = {'-', '-', '-'};
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 3]);

    %% Plot all ERPs together
    subplot(numel(cond_lab),1,1); hold on;
    plot_info.title  = strjoin(cond_lab,',');
    plot_info.ax     = gca;
    
    % Compute means and variance
    means = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
    sem   = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = roi_erp{cond_ix}.avg(ch_ix,:);
        sem(cond_ix,:) = squeeze(sqrt(roi_erp{cond_ix}.var(ch_ix,:))./sqrt(numel(roi_erp{cond_ix}.cfg.previous.trials)))';
    end
    
    fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
    
    %% Plot Stat Comparisons
    for st_ix = 1:numel(cond_comp)
        cond_ixs = cond_comp{st_ix};
        subplot(numel(cond_lab),1,st_ix+1); hold on;
        plot_info.title  = strjoin(cond_lab(cond_ixs),'-');
        plot_info.ax     = gca;
        % Condition plotting params
        cond_info.name       = cond_lab(cond_ixs);
        cond_info.style      = repmat({'-'},size(cond_ixs));
        cond_info.color      = cond_colors(cond_ixs);
        cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_ixs)]);
        
        % Compute means and variance
        means = NaN([numel(cond_ixs) size(roi_erp{1}.avg,2)]);
        sem   = NaN([numel(cond_ixs) size(roi_erp{1}.avg,2)]);
        for c_ix = 1:numel(cond_ixs)
            means(c_ix,:) = roi_erp{cond_ixs(c_ix)}.avg(ch_ix,:);
            sem(c_ix,:)   = squeeze(sqrt(roi_erp{cond_ixs(c_ix)}.var(ch_ix,:))./sqrt(numel(roi_erp{cond_ixs(c_ix)}.cfg.previous.trials)))';
        end
        
        %Find significant time periods
        if sum(stat{st_ix}.mask(ch_ix,:))>0
            sig_chunks = fn_find_chunks(stat{st_ix}.mask(ch_ix,:));
            sig_chunks(stat{st_ix}.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
            % If stat and roi_erp aren't on same time axis, adjust sig_chunk indices
            if (size(stat{st_ix}.time,2)~=size(roi_erp{1}.time,2)) || (sum(stat{st_ix}.time==roi_erp{1}.time)~=numel(stat{st_ix}.time))
                for chunk_ix = 1:size(sig_chunks,1)
                    sig_chunks(chunk_ix,1) = find(roi_erp{1}.time==stat{st_ix}.time(sig_chunks(chunk_ix,1)));
                    sig_chunks(chunk_ix,2) = find(roi_erp{1}.time==stat{st_ix}.time(sig_chunks(chunk_ix,2)));
                end
            end
            fprintf('%s -- %i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                stat{st_ix}.label{ch_ix},size(sig_chunks,1));
            fn_plot_ts_errbr_sig(plot_info,means,sem,sig_chunks,event_info,cond_info);
        else
            fprintf('%s -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',stat{st_ix}.label{ch_ix});
            fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
        end
    end
    
    %% Save stats variable
    stats_fname = [SBJ_vars.dirs.proc SBJ '_' an_id '.mat'];
    save(stats_fname, '-v7.3', 'design', 'roi_erp', 'stat');
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_filename]);
    end
end
end

