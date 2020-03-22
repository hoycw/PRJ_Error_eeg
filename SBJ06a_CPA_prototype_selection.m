function SBJ06a_CPA_prototype_selection(SBJ, proc_id, cpa_id)
%% Candidate-Prototype Analysis (CPA): Prototype Selection
%   Selects top IC based on spatial (elec_list) and temporal (time_win)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
cpa_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' cpa_id '_vars.m'];
eval(cpa_vars_cmd);

%% Load the data
%loaded data from after SBJ02a --> already cleaned and trial segmented
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat']); %chose 02a - ica before rejection!
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

[cond_lab, ~, cond_colors, cond_styles, ~] = fn_condition_label_styles('Odd'); % maybe change this so not hardcoded
cond_idx = fn_condition_index(cond_lab, bhv);

% Create contrast: (Unexpected - Expected) for each outcome
[diff_lab, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles(cpa.diff_id);
if numel(diff_lab)>1; error('Too many condition contrasts!'); end

%% Exclude bad trials
% Find trials that overlap with bad_epochs from raw visual inspection
load([SBJ_vars.dirs.preproc SBJ '_preproc_eeg_full_ft.mat'],'data');
load([SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat']);
if ~isempty(bad_epochs)
    bad_raw_trials = fn_find_trials_overlap_epochs(bad_epochs,1:size(data.trial{1},2),...
        event_onsets,proc.trial_lim_s*data.fsample);    
else
    bad_raw_trials = [];
end

% Identify training and bad behavioral trials
training_ix = find(bhv.blk==0);
rt_low_ix   = find(bhv.rt <= proc.rt_bounds(1) & bhv.rt>0);
rt_high_ix  = find(bhv.rt >= proc.rt_bounds(2));
%exclude_trials = unique(vertcat(bad_raw_trials, training_ix));
exclude_trials = unique(vertcat(bad_raw_trials, training_ix, rt_low_ix, rt_high_ix));

% Exclude bad trials
cfgs = [];  
cfgs.trials = setdiff([1:numel(trials.trial)], exclude_trials');
trials = ft_selectdata(cfgs, trials);
eog_trials = ft_selectdata(cfgs, eog_trials);

%% Temporal Criteria: Condition Stats
% Select stats epoch
cfg = [];
% NOT THIS SIMPLE: cfg.trials  = setdiff(1:numel(ica.trial),SBJ_vars.trial_reject_ix);
%   need to toss orig_exclude (bad_epochs, RTs), then trial_reject_ix
cfg.latency = cpa.time_win;
st_ica = ft_selectdata(cfg,ica);

st_trials = cell([2 1]);
for d_ix = 1:2
    cond_trial_ix = find(cond_idx==diff_pairs{1}(d_ix));
    st_trials{d_ix} = nan([numel(st_ica.label) numel(cond_trial_ix) numel(st_ica.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        st_trials{d_ix}(:,t_ix,:) = st_ica.trial{cond_trial_ix(t_ix)};
    end
end

% diff_waves = zeros(numel(st_data.label), numel(diff_lab), numel(st_data.time{1}));
sig_wins = nan([numel(ica.label) numel(st_ica.time{1})]);
p_vals   = nan([numel(ica.label) numel(st_ica.time{1})]);
sig_perc = zeros(size(ica.label));
sig_len  = zeros(size(ica.label));
for comp_ix = 1:numel(st_ica.label)
    % Compute stats between conditions in contrast within time windows
    for time_ix = 1:numel(st_ica.time{1})
        [sig_wins(comp_ix, time_ix), p_vals(comp_ix, time_ix)] = ...
            ttest2(st_trials{1}(comp_ix, :, time_ix), ...
                   st_trials{2}(comp_ix, :, time_ix));
    end
    if cpa.alpha~=0.05; error('re-compute significance with non 0.05 threshold'); end
    
    % Compute summary metrics (% sig, min_sig_len)
    sig_perc(comp_ix) = sum(sig_wins(comp_ix,:)) / size(sig_wins,2);
    [sig_lims] = fn_find_chunks(sig_wins(comp_ix,:));
    sig_lims(squeeze(sig_wins(comp_ix,sig_lims(:,1)))==0,:) = [];
    if ~isempty(sig_lims)
        sig_len(comp_ix) = max(diff(sig_lims,1,2)+1);
    end
    
%     %Difference Wave (for debugging)
%     diff_waves(comp_ix, diff_ix,:) = plot_means(diff_ix,:); % whole time period
end

% Select components with consecutive significance
time_ic_idx = sig_len./st_ica.fsample >= cpa.min_sig_len;
fprintf('%s: %d / %d components meet temporal criteria!\n',SBJ,sum(time_ic_idx),numel(ica.label));

%% Spatial Criteria: Topo Matching
% Average ERP into topo map
%erp_topos = nan([numel(ica.label) numel(ica.topolabel)]);
% !!!Implement grubbs' test???
%   can use isoutlier(data,'method','grubbs');
space_ic_idx = false(size(ica.label));
if strcmp(cpa.elec_method,'peak')
    top_elecs = cell([numel(ica.label) cpa.n_max_elec]);
    for comp_ix = 1:numel(ica.label)
        [~, top_ix] = maxk(abs(ica.topo(:,comp_ix)), cpa.n_max_elec);
        top_elecs(comp_ix,:) = ica.topolabel(top_ix);
        if numel(intersect(top_elecs(comp_ix,:),cpa.elec_list)) >= cpa.min_elec_match
            space_ic_idx(comp_ix) = true;
        end
    end
end
fprintf('%s: %d / %d components meet spatial criteria!\n',SBJ,sum(space_ic_idx),numel(ica.label));

%% Combine Criteria
var_ic_idx = true(size(ica.label));
if isfield(cpa,'ic_rank_max')
    % Select components ranked highest when ordered by variance
    var_ic_idx(cpa.ic_rank_max+1:end) = false;
end
good_ic_idx = true(size(ica.label));
good_ic_idx(SBJ_vars.ica_reject) = false;
fprintf('%s: %d / %d components were kept in cleaning!\n',SBJ,sum(good_ic_idx),numel(good_ic_idx));

final_ics = find(all([time_ic_idx, space_ic_idx, var_ic_idx, good_ic_idx],2));
fprintf('%s: %d / %d components selected!\n',SBJ,numel(final_ics),numel(ica.label));
if sum(final_ics)<1
    error('No ICs found that match all criteria!');
end

%% Save Data
clean_data_fname = [SBJ_vars.dirs.proc SBJ '_' cpa_id '_' proc_id '_prototype.mat'];
save(clean_data_fname, '-v7.3', 'final_ics');

end
