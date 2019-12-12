function SBJ06a_CPA_TT(SBJ, proc_id, plt_id, electrodes, time_win)

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/', ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load the data
%loaded data from after SBJ02a --> already cleaned and trial segmented
load([SBJ_vars.dirs.preproc SBJ '_preproc_eeg_full_ft.mat']);
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat']); %chose 02a - ica before rejection!
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

clean_ica = ica;
%% trim data to plotting time -- is this required?
cfg = [];
%cfg.latency = plt.plt_lim; i dont think this is necessary? (already -.2 to
%1.3)
data = ft_selectdata(cfg,clean_ica);
time_vec = data.time{1};
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles('DifOutS'); % maybe change this so not hardcoded
cond_idx = fn_condition_index(cond_lab, bhv);
% Get trials for plotting
trials = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(data.label) numel(cond_trial_ix) numel(data.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,t_ix,:) = data.trial{cond_trial_ix(t_ix)};
    end
end
[~, min_ix] = min(abs(clean_ica.time{1,1}(:) - time_win(1)));
[~, max_ix] = min(abs(clean_ica.time{1,1}(:) - time_win(2)));
for comp_ix = 1:numel(data.label)
    %% Compute plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) max_ix-min_ix+1]);
    sems  = NaN([numel(cond_lab) max_ix-min_ix+1]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(comp_ix,:,min_ix:max_ix),2))';
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(comp_ix,:,min_ix:max_ix),[],2))./sqrt(size(trials{cond_ix},2))';
        means_all(comp_ix, cond_ix, :) = means(cond_ix,:);
        sems_all(comp_ix, cond_ix, :) = sems(cond_ix,:);
    end
%avg_time_win(ch_ix) = mean(abs(diff_waves(ch_ix, min_ix: max_ix))); %not necessary anymore -- calculates average amplitude of dfifference wave in given time window
end

%% compute ERP for each Electrode

%% PLOT AND SAVE
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/Dec_11_Graphs/' SBJ '/'] ;
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

sig_chunks = cell(numel(data.label), 1);
for comp_ix = 1: numel(data.label)
    a = '';
    for elec_ix = 1:numel(electrodes)
        a = [a electrodes{elec_ix} '_'];
    end
    fig_name = [SBJ 'Component0' num2str(comp_ix) 'Electrodes_' a '_TT_'];
    figure('Name', fig_name, 'units','normalized',...
        'outerposition',[0 0 0.5 0.8], 'Visible', 'off');
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+1 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(time_vec(min_ix:max_ix), means_all(comp_ix, cond_ix, :), sems_all(comp_ix, cond_ix,:),...
            {'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},plt.errbar_alpha);
        hold on
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    ax.YLabel.String = 'uV';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    leg_lab = [cond_lab 'F' cond_lab(cond_ix)];
    %if plt.legend
       % legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
   % end
     fig_fname = [fig_dir fig_name '.png'];
     fprintf('Saving %s\n',fig_fname);
     % Ensure vector graphics if saving
     saveas(gcf,fig_fname);
end    

end