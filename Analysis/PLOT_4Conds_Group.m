function PLOT_4Conds_Group(SBJs, plt_id, an_id, fig_vis, save_fig, fig_ftype)
%Purpose: This function loads the data from PLOT_4Conds_Individ and computes a total average and plots it for all 4 conditons (easy/miss, easy/hit, hard/hit, hard/miss).
%SBJS = cell of strings of subjects you want to plot {'EEG01', 'EEG02',
%etc}
%proc_id = 
%plt_id = 'ts_F15to28_evnts_sigPatch'
%an_id = 'ERP_Cz_F_trl15t28_flt05t20_stat06'
%fig_vis = whether or not the plot should pop up ('on' or 'off')
%save_fig = save  the figure or not (1 or 0)
%fig_ftype = file type the figure should be saved as (i.e. 'png')
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
cond_lab_hit = [1 0];
cond_lab = {'easy', 'hard'};
cond_colors = {[0.6350, 0.0780, 0.1840], [0.3010, 0.7450, 0.9330], 	[0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560]};
stats_dir = [root_dir 'PRJ_Error_eeg/results/Stats/4_Conditions/' an_id '/'];
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/4Conds/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
for s = 1:length(SBJs)
    clean_stats_fname{s} = [stats_dir an_id SBJs{s} '_4Conds.mat'];
    load(clean_stats_fname{s});
end

% load all SBJ roi_erp
roi_erp = cell(size(SBJs));
for s = 1:length(SBJs)
    tmp = load(clean_stats_fname{s}); roi_erp{s} = tmp.roi_erp;
end
num_conds = (size(cond_lab, 2)+size(cond_lab_hit, 2));
gavg = cell(size(num_conds));
averages = cell(1, numel(SBJs));
for cond_ix = 1:(numel(cond_lab)+numel(cond_lab_hit))
    for index = 1: numel(SBJs)
        averages{index}= roi_erp{1,index}{cond_ix};
    end
    gavg{cond_ix} = ft_timelockgrandaverage(cfg_avg, averages{1:numel(SBJs)});
end

%%
for ch_ix = 1:numel(an.ROI)
    fig_name = ['4Conds_Avg' an.ROI{ch_ix}];
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
    [~,event_info.time]  = min(abs(roi_erp{1}{1}.time-0.18));
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
    plot_info.title  = 'Group 4 Conditions';
    plot_info.ax     = gca;
%% calculate statistics
        % Compute means and variance
        means = NaN([numel(cond_lab)+numel(cond_lab_hit) size(gavg{1}.avg,2)]);
        sem   = NaN([numel(cond_lab)+numel(cond_lab_hit) size(gavg{1}.avg,2)]);
        for c_ix = 1:numel(cond_lab)+numel(cond_lab_hit)
            means(c_ix,:) = gavg{c_ix}.avg(ch_ix,:);
            sem(c_ix,:) = squeeze(sqrt(gavg{c_ix}.var(ch_ix,:))./sqrt(numel(SBJs))); %What do I divide by? Whats the population size?
        end
        
    
    fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
    
   
    
    %% Save stats variable
    stats_fname = [stats_dir an_id 'Group_4Conds.mat'];
    save(stats_fname, '-v7.3', 'roi_erp', 'gavg');
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_filename]);
    end
    
end
end
