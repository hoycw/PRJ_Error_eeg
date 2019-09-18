function PLOT_Indiv_4Conds(SBJ, proc_id, plt_id, an_id, fig_vis, save_fig, fig_ftype)

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

%% Compute ERPs
% Select Channel(s)
cfgs = [];
cfgs.channel = an.ROI;
roi = ft_selectdata(cfgs, clean_trials);

cond_lab = {'easy', 'hard'};
cond_lab_hit = [1 0];
cond_colors = {[0.6350, 0.0780, 0.1840], [0.3010, 0.7450, 0.9330], 	[0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560]}

num_conds = (size(cond_lab, 2)+size(cond_lab_hit, 2));
roi_erp  = cell(1,num_conds);
n_trials = zeros(1, num_conds);
cfgavg = [];
cfgavg.keeptrials = 'no';
cfgavg.paramater = 'avg';

for cond_ix = 1:numel(cond_lab)
    for cond_ix_hit = 1:numel(cond_lab_hit);
        trials_hit = find(bhv.hit == cond_lab_hit(cond_ix_hit));
        trials_cond = find(strcmp(bhv.cond, cond_lab{cond_ix}));
        cfgavg.trials = intersect(trials_hit, trials_cond);
        index = (cond_ix*2 - 1)+(cond_ix_hit - 1);
        roi_erp{index} = ft_timelockanalysis(cfgavg,roi);
        n_trials(index) = size(cfgavg.trials, 1);
    end
end

%{
for st_ix = 1
    cond_ixs = [1 2 3 4];
    
    % Create make design matrix and stats
    design = zeros(2, sum(n_trials(cond_ixs(:))));
    for c_ix = 1:4
        if c_ix==1
            design(1,1:n_trials(cond_ixs(c_ix))) = cond_ixs(c_ix);                               % Conditions (Independent Variable)
            design(2,1:n_trials(cond_ixs(c_ix))) = 1:n_trials(cond_ixs(c_ix));
            previous = n_trials(cond_ixs(c_ix));
            % Trial Numbers
        else
            design(1,sum(n_trials(cond_ixs(1:c_ix-1)))+1:(previous+n_trials(cond_ixs(c_ix)))) = cond_ixs(c_ix); % Conditions (Independent Variable)
            design(2,sum(n_trials(cond_ixs(1:c_ix-1)))+1:(previous+n_trials(cond_ixs(c_ix)))) = 1:n_trials(cond_ixs(c_ix));
            previous = previous + n_trials(cond_ixs(c_ix));
        end
    end
    
    % Compute stats between ERPs
    cfg_stat.design           = design;
    [stat{st_ix}] = ft_timelockstatistics(cfg_stat, roi_erp{cond_ixs(1)}, roi_erp{cond_ixs(2)});
end
%}

if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' SBJ '/4_Conditions/' an_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
end
for ch_ix = 1:numel(an.ROI)
    fig_name = [SBJ '_4Conds_' an.ROI{ch_ix}];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible', fig_vis);   %this size is for single plots
    
    % General plotting params
    plot_info.fig        = f;
    plot_info.sig_color  = plt_vars.sig_color;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.x_step     = plt_vars.x_step_sz*roi.fsample;
    plot_info.x_lab      = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    plot_info.y_lab      = 'uV';
    plot_info.legend     = plt_vars.legend;
    plot_info.legend_loc = plt_vars.legend_loc;
    % Event plotting params
    event_info.name      = {an.event_type};
    [~,event_info.time]  = min(abs(roi_erp{1}.time-1.8));
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = {plt_vars.evnt_color};
    event_info.style     = {plt_vars.evnt_style};
    % Condition plotting params
    cond_info.name       = {'Correct/Easy', 'Error/Easy', 'Correct/Hard', 'Error/Hard'};
    cond_info.style      = {'-', '-', '-', '-'};
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 4]);

    %% Plot all ERPs together
    hold on;
    plot_info.title  = [SBJ ' 4 Conditions'];
    plot_info.ax     = gca;
    
    % Compute means and variance
    means = NaN([numel(cond_lab)+numel(cond_lab_hit) size(roi_erp{1}.avg,2)]);
    sem   = NaN([numel(cond_lab)+numel(cond_lab_hit) size(roi_erp{1}.avg,2)]);
    for cond_ix = 1:numel(cond_lab)+numel(cond_lab_hit);
        means(cond_ix,:) = roi_erp{cond_ix}.avg(ch_ix,:);
        sem(cond_ix,:) = squeeze(sqrt(roi_erp{cond_ix}.var(ch_ix,:))./sqrt(numel(roi_erp{cond_ix}.cfg.previous.trials)))';
    end
    
    fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
    if save_fig
        stats_fname = [SBJ_vars.dirs.proc SBJ '_' an_id 'stats_4cond.mat'];
        save(stats_fname, '-v7.3', 'roi_erp');
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_filename]);
    end
end