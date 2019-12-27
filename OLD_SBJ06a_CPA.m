function OLD_SBJ06a_CPA(SBJ, proc_id, plt_id, electrodes, time_win)

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
cfg.latency = plt.plt_lim;
data = ft_selectdata(cfg,clean_ica);

[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles('Odd'); % maybe change this so not hardcoded
cond_idx = fn_condition_index(cond_lab, bhv);
% Create contrast: (Unexpected - Expected) for each outcome
[diff_lab, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles('DiffOutStdTar');
% Get trials for plotting
trials = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(data.topolabel) numel(cond_trial_ix) numel(data.time{1})]);
    for trial_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,trial_ix,:) = data.trial{cond_trial_ix(trial_ix)};
    end
end

%% compute ERP
cfg = [];
cfg.channel = 'all';
erps = ft_timelockanalysis(cfg, data);
components_ica_reject = unique([SBJ_vars.ica_reject, heog_ics, veog_ics]);
components = 1:numel(clean_ica.topo(1,:));
components_keep = setdiff(components, components_ica_reject);
for comp_ix = 1:numel(components_keep)
    avg_topo(comp_ix) = mean(clean_ica.topo(:,comp_ix)); % do i want average for all topos? average for a component?
end
for elec_ix = 1:numel(electrodes)
    temp_elec = electrodes{elec_ix}
    for topo_elec_ix = 1:numel(clean_ica.topolabel)
        temp_ica = clean_ica.topolabel{topo_elec_ix};
        if (strcmp(temp_elec, temp_ica))
            index_electrodes(elec_ix) = topo_elec_ix;
        end
    end
end
topo_components = [];
for elec_ix = 1:numel(index_electrodes)
    [temp, ~] = find(max == index_electrodes(elec_ix));
    topo_components = union(topo_components, temp);
end
topo_components = topo_components';

diff_waves = zeros(numel(data.topolabel), numel(data.time{1}));
[~, min_ix] = min(abs(clean_ica.time{1,1}(:) - time_win(1)));
[~, max_ix] = min(abs(clean_ica.time{1,1}(:) - time_win(2)));
for ch_ix = 1:numel(data.topolabel)
    %% Compute plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(data.time{1})]);
    sems  = NaN([numel(cond_lab) numel(data.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(ch_ix,:,:),2));
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(ch_ix,:,:),[],2))./sqrt(size(trials{cond_ix},2))';
    end
    
    plot_means = NaN([numel(diff_lab) numel(data.time{1})]);
    for diff_ix = 1:numel(diff_lab)
        if iscell(diff_pairs{diff_ix})
            contrast_pair = diff_pairs{diff_ix}{2};
            contrast_mean = mean(means(contrast_pair{:},:),1);
            plot_means(diff_ix,:) = means(cond_ix,:)-contrast_mean;
        elseif any(diff_pairs{diff_ix}==0)
            cond_ix = diff_pairs{diff_ix}(diff_pairs{diff_ix}~=0);
            plot_means(diff_ix,:) = means(cond_ix,:);
        else
            [sig(ch_ix), p(ch_ix)] = ttest(means(diff_pairs{diff_ix}(1),min_ix:max_ix), means(diff_pairs{diff_ix}(2),min_ix:max_ix)); % tests the distribution of cond1 vs. cond2
            plot_means(diff_ix,:) = means(diff_pairs{diff_ix}(1),:)-means(diff_pairs{diff_ix}(2),:);
        end
    end
diff_waves(ch_ix,:) = plot_means(diff_ix,:);
avg_time_win(ch_ix) = mean(abs(diff_waves(ch_ix, min_ix: max_ix))); %not necessary anymore -- calculates average amplitude of dfifference wave in given time window
end
ica_reject = intersect(unique([SBJ_vars.ica_reject, heog_ics, veog_ics]), topo_components);
components_final_no_ica = setdiff(topo_components, ica_reject);
[~, components_temp] = mink(p(topo_components), 3);
components = topo_components(components_temp);
%[~,erp_components] = find(sig == 1);

%components = intersect(topo_components, erp_components);
%ica_reject = intersect(unique([SBJ_vars.ica_reject, heog_ics, veog_ics]), components);
%components_final_no_ica = setdiff(components, ica_reject);

%% Save Data
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '_06a.mat'];
save(clean_data_fname, '-v7.3', 'components_final_no_ica', 'components');
end
