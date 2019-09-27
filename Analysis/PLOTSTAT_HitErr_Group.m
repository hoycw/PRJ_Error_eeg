function PLOTSTAT_HitErr_Group(SBJs, plt_id, an_id, fig_vis, save_fig, fig_ftype)
%Purpose:This function loads the subjects data retrieved from running PLOTSTAT_HitErr_Invid.  Then it plots the ERPs at the group level for the hit vs miss conditions, with error bars and significant clusters.
%Inputs
%SBJ = string (ie) 'EEG01'
%proc_id = 'eeg_full_ft'
%plt_id = 'ts_F15to28_evnts_sigPatch'
%an_id = 'ERP_Cz_F_trl15t28_flt05t20_stat06'
%fig_vis = whether or not the plot should pop up ('on' or 'off')
%save_fig = 0 or 1
%fig_ftype = 'png' (a string)

%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);
addpath(ft_dir);
ft_defaults

%% Load preprocessed data
for s = 1:length(SBJs)
    SBJ_vars_cmd{s} = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd{s});
    SBJ_vars_all{s} = SBJ_vars;
end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cond_lab = [1 0];
cond_colors = {[0 0 1],[0 1 0],[1 0 0]};
stats_dir = [root_dir 'PRJ_Error_eeg/results/Stats/HitErr/' an_id '/'];
fig_dir = [root_dir 'PRJ_Error_eeg/results/Stats/HitErr/' an_id '/'];
if ~exist(fig_dir,'dir')
   mkdir(fig_dir);
end
for s = 1:length(SBJs)
    clean_stats_fname{s} = [stats_dir an_id SBJs{s} '_HitErr.mat'];
    load(clean_stats_fname{s});
end

% load all SBJ roi_erp
roi_erp = cell(size(SBJs));
for s = 1:length(SBJs)
    tmp = load(clean_stats_fname{s}); roi_erp{s} = tmp.roi_erp;
end

gavg = cell(size(cond_lab));
averages = cell(numel(SBJs));
for cond_ix = 1:numel(cond_lab)
    for index = 1: numel(SBJs)
        averages{index}= roi_erp{1, 1}{cond_ix};
    end
    gavg{cond_ix} = ft_timelockgrandaverage(cfg_avg, averages{1:numel(SBJs)});
end

stat   = cell(size(1));
design = cell(size(1));
n_subjs = numel(SBJs);
for st_ix = 1
    cond_ixs = [1 2];
    
    % Create make design matrix and stats
    design{st_ix} = zeros(2, n_subjs);
    for c_ix = 1:2
        if c_ix==1
            design{st_ix}(1,1:n_subjs) = cond_ixs(c_ix);                               % Conditions (Independent Variable)
            design{st_ix}(2,1:n_subjs) = 1:n_subjs;                   % Trial Numbers
        else
            design{st_ix}(1,n_subjs+1:n_subjs*c_ix) = cond_ixs(c_ix); % Conditions (Independent Variable)
            design{st_ix}(2,(n_subjs+1):n_subjs*c_ix)= 1:n_subjs;
        end
    end
    
    % Compute stats between ERPs
    cfg_stat.design           = design{st_ix};
    cfg_stat.parameter        = 'avg';
    [stat{st_ix}] = ft_timelockstatistics(cfg_stat, gavg{cond_ixs(1)}, gavg{cond_ixs(2)});
end

%%
for ch_ix = 1:numel(an.ROI)
    fig_name = ['HitErr_Avg' an.ROI{ch_ix}];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible', fig_vis);   %this size is for single plots
    %%
    % General plotting params
    plot_info.fig        = f;
    plot_info.sig_color  = plt_vars.sig_color;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.x_step = plt_vars.x_step_sz*plt_vars.fsample;
    plot_info.x_lab  = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    plot_info.y_lab  = 'uV';
    plot_info.legend = plt_vars.legend;
    plot_info.legend_loc = plt_vars.legend_loc;
    % Event plotting params
    event_info.name      = {an.event_type};
    %%
    [~,event_info.time]  = min(abs(roi_erp{1}{1}.time-0));
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = {plt_vars.evnt_color};
    event_info.style     = {plt_vars.evnt_style};
    % Condition plotting params
    cond_info.name       = {'correct', 'error'};
    cond_info.style      = {'-', '-', '-'};
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 3]);

    %% Plot all ERPs together
    hold on;
    plot_info.title  = 'Easy vs. Hard Conditions';
    plot_info.ax     = gca;
%% calculate statistics
        % Compute means and variance
        means = NaN([numel(cond_lab) size(gavg{1}.avg,2)]);
        sem   = NaN([numel(cond_lab) size(gavg{1}.avg,2)]);
        temp = cell(size(numel(SBJs)));
        store = zeros(numel(SBJs),length(gavg{1}.avg));
        for c_ix = 1:numel(cond_lab)
            means(c_ix,:) = gavg{c_ix}.avg(ch_ix,:);
            for s = 1:numel(SBJs)
                temp{s} = roi_erp{s}{c_ix};
                for n = 1: length(gavg{1}.avg)
                    store(s,n) = temp{s}.avg(n);
                end
            end
            for t_ix = 1: length(gavg{1}.avg)
                sem(c_ix,t_ix) = std(store(:,t_ix))/ sqrt(numel(SBJs));
            end
        end
        
    
    fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
    
    %% Plot Stat Comparisons
    for st_ix = 1
        cond_ixs =[1 2];
        subplot(numel(cond_lab),1,st_ix); hold on;
        plot_info.title  = 'Easy vs. Hard Conditions';
        plot_info.ax     = gca;
        % Condition plotting params
        cond_info.name       = {'correct', 'error'};
        cond_info.style      = repmat({'-'},size(cond_ixs));
        cond_info.color      = cond_colors(cond_ixs);
        cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_ixs)]);

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
    stats_fname = [fig_dir an_id 'Group_HitErr.mat'];
    save(stats_fname, '-v7.3', 'design', 'roi_erp', 'stat', 'gavg');
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_filename]);
    end
end
end

