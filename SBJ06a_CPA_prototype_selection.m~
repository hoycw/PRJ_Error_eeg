function SBJ06a_CPA_prototype_selection(SBJ, proc_id, cpa_id,save_fig,varargin)
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

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
cpa_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' cpa_id '_vars.m'];
eval(cpa_vars_cmd);

%% Load the data
%loaded data from after SBJ02a --> already cleaned and trial segmented
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat'],'ica'); %chose 02a - ica before rejection!
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'],'bhv');
load([SBJ_vars.dirs.events SBJ '_' proc_id '_02a_orig_exclude_trial_ix.mat']);

[cond_lab, ~, cond_colors, cond_styles, ~] = fn_condition_label_styles('Odd'); % maybe change this so not hardcoded
cond_idx = fn_condition_index(cond_lab, bhv);

% Create contrast: (Unexpected - Expected) for each outcome
[diff_lab, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles(cpa.diff_id);
if numel(diff_lab)>1; error('Too many condition contrasts!'); end

%% Exclude bad trials
cfgs = [];  
cfgs.trials = setdiff([1:numel(ica.trial)], SBJ_vars.trial_reject_ix_oddball);
clean_ica = ft_selectdata(cfgs, ica);

%% Temporal Criteria: Condition Stats
% Select stats epoch
cfg = [];
cfg.latency = cpa.time_win;
st_ica = ft_selectdata(cfg,clean_ica);

st_trials = cell([2 1]);
for d_ix = 1:2
    cond_trial_ix = find(cond_idx==diff_pairs{1}(d_ix));
    st_trials{d_ix} = nan([numel(st_ica.label) numel(cond_trial_ix) numel(st_ica.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        st_trials{d_ix}(:,t_ix,:) = st_ica.trial{cond_trial_ix(t_ix)};
    end
end

% diff_waves = zeros(numel(st_data.label), numel(diff_lab), numel(st_data.time{1}));
sig_wins = nan([numel(st_ica.label) numel(st_ica.time{1})]);
p_vals   = nan([numel(st_ica.label) numel(st_ica.time{1})]);
sig_perc = zeros(size(st_ica.label));
sig_len  = zeros(size(st_ica.label));
sig_len_chunk  = zeros(size(st_ica.label));
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
        sig_len(comp_ix) = sum(diff(sig_lims,1,2)+1);
        sig_len_chunk(comp_ix) = max(diff(sig_lims,1,2)+1);
    end
    
%     %Difference Wave (for debugging)
%     diff_waves(comp_ix, diff_ix,:) = plot_means(diff_ix,:); % whole time period
end

% Convert to time
sig_len = sig_len./st_ica.fsample;
sig_len_chunk = sig_len_chunk./st_ica.fsample;

% Select components with consecutive significance
time_ic_idx = sig_len >= cpa.min_sig_len;
fprintf('%s: %d / %d components meet temporal criteria!\n',SBJ,sum(time_ic_idx),numel(st_ica.label));

%% Spatial Criteria: Topo Matching
% Compute topo outliers
%   isoutlier(data,'method','grubbs');

% Compute top electrodes
top_elecs = cell([numel(st_ica.label) cpa.n_max_elec]);
n_elec_match = zeros(size(st_ica.label));
for comp_ix = 1:numel(st_ica.label)
    [~, top_ix] = maxk(abs(clean_ica.topo(:,comp_ix)), cpa.n_max_elec);
    top_elecs(comp_ix,:) = clean_ica.topolabel(top_ix);
    n_elec_match(comp_ix) = numel(intersect(top_elecs(comp_ix,:),cpa.elec_list));
end

% Compute topo correlations
% Load and average group topo
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_final.mat'],'clean_trials'); %chose 02a - ica before rejection!
load([root_dir 'PRJ_Error_eeg/data/GRP/' cpa.topo_SBJ_id '_Odd_' cpa.topo_an_id '.mat']);
cfg_avg = [];
cfg_avg.channel = clean_trials.label;
cfg_avg.latency = cpa.topo_lim;
cfg_avg.avgovertime = 'yes';
topo = ft_selectdata(cfg_avg,er_grp{strcmp(cond_lab,cpa.topo_cond)});
topo_match = true(size(topo.label));
if numel(clean_ica.topolabel)~=numel(topo.label) || ~all(strcmp(clean_ica.topolabel,topo.label))
    topo_order = zeros(size(topo.label));
    for lab_ix = 1:numel(topo.label)
        if any(strcmp(clean_ica.topolabel,topo.label{lab_ix}))
            topo_order(lab_ix) = find(strcmp(clean_ica.topolabel,topo.label{lab_ix}));
        else
            topo_match(lab_ix) = false;
        end
    end
    topo_order(topo_order==0) = [];
else
    topo_order = 1:numel(topo.label);
end

% Correlate ICA and grand average topos
topo_corrs = nan(size(st_ica.label));
topo_pvals = nan(size(st_ica.label));
for comp_ix = 1:numel(st_ica.label)
    [tmp_corr, tmp_pval] = corrcoef(clean_ica.topo(topo_order,comp_ix), topo.avg(topo_match));
    topo_corrs(comp_ix) = tmp_corr(1,2);
    topo_pvals(comp_ix) = tmp_pval(1,2);
end

% Combine spatial criteria
if any(strcmp(cpa.elec_method,'peak'))
	peak_idx = n_elec_match >= cpa.min_elec_match;
else
    peak_idx = true;
end
if any(strcmp(cpa.elec_method,'topo_corr'))
	topo_corr_idx = topo_pvals <= cpa.topo_pval;
else
    topo_corr_idx = true;
end
space_ic_idx = all([peak_idx, topo_corr_idx],2);
fprintf('%s: %d / %d components meet peak criterion!\n',SBJ,sum(peak_idx),numel(clean_ica.label));
fprintf('%s: %d / %d components meet topo corr criterion!\n',SBJ,sum(topo_corr_idx),numel(clean_ica.label));
fprintf('%s: %d / %d components meet spatial criteria!\n',SBJ,sum(space_ic_idx),numel(clean_ica.label));

%% Combine Criteria
var_ic_idx = true(size(clean_ica.label));
if isfield(cpa,'ic_rank_max')
    % Select components ranked highest when ordered by variance
    var_ic_idx(cpa.ic_rank_max+1:end) = false;
end
good_ic_idx = true(size(clean_ica.label));
good_ic_idx(SBJ_vars.ica_reject) = false;
fprintf('%s: %d / %d components were kept in cleaning!\n',SBJ,sum(good_ic_idx),numel(good_ic_idx));

final_ics = find(all([time_ic_idx, space_ic_idx, var_ic_idx, good_ic_idx],2));
fprintf('%s: %d / %d components selected!\n',SBJ,numel(final_ics),numel(clean_ica.label));
if sum(final_ics)<1
    warning('No ICs found that match all criteria!');
end

%% QA plots
% Plot individual criteria
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/prototype/' cpa_id '/QA_all/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

fig_name = [SBJ '_proto_' cpa_id '_QA_all'];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
ns_sz = 25;
sig_sz = 75;

% Temporal Cirterion
subplot(3,1,1); hold on;
ylim([0 diff(cpa.time_win)]);
scatter(1:numel(clean_ica.label),sig_len,ns_sz,'k');
scatter(find(time_ic_idx),sig_len(time_ic_idx),sig_sz,'k','filled');
scatter(final_ics,sig_len(final_ics),sig_sz,'r','filled');
line(xlim, [cpa.min_sig_len cpa.min_sig_len],'Color','r');
line([cpa.ic_rank_max cpa.ic_rank_max],ylim,'Color','r');
xlabel('IC #');
ylabel('min sig len (s)');
title([num2str(sum(time_ic_idx)) '/' num2str(numel(final_ics)) '; Window: ' ...
    num2str(cpa.time_win(1)) '-' num2str(cpa.time_win(2)) ' (' num2str(diff(cpa.time_win)) ' s)']);
set(gca,'FontSize',14);

% Spatial Peak Elecs Criterion
subplot(3,1,2); hold on;
ylim([0 5]);
scatter(1:numel(clean_ica.label),n_elec_match,ns_sz,'k');
scatter(find(peak_idx),n_elec_match(peak_idx),ns_sz,'k','filled');
scatter(find(space_ic_idx),n_elec_match(space_ic_idx),sig_sz,'k','filled');
scatter(final_ics,n_elec_match(final_ics),sig_sz,'r','filled');
line(xlim, [cpa.min_elec_match cpa.min_elec_match],'Color','r');
line([cpa.ic_rank_max cpa.ic_rank_max],ylim,'Color','r');
ylabel('# elec match');
title([num2str(sum(peak_idx)) '/' num2str(numel(final_ics)) '; elecs: ' strjoin(cpa.elec_list)]);
xlabel('IC #');
set(gca,'FontSize',14);

% Spatial Topo Corr Criterion
subplot(3,1,3); hold on;
ylim([-1 1]);
scatter(1:numel(clean_ica.label),topo_corrs,ns_sz,'k');
scatter(find(topo_corr_idx),topo_corrs(topo_corr_idx),ns_sz,'k','filled');
scatter(find(space_ic_idx),topo_corrs(space_ic_idx),sig_sz,'k','filled');
scatter(final_ics,topo_corrs(final_ics),sig_sz,'r','filled');
% line(xlim, [cpa.min_elec_match cpa.min_elec_match],'Color','r');
line([cpa.ic_rank_max cpa.ic_rank_max],ylim,'Color','r');
ylabel('topo corr');
title([num2str(sum(topo_corr_idx)) '/' num2str(numel(final_ics)) '; ' ...
    cpa.topo_cond ' topo (' num2str(cpa.topo_lim(1)) '-' num2str(cpa.topo_lim(2)) ' s)']);
xlabel('IC #');
set(gca,'FontSize',14);

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot combined spatial criteria (Peak vs. topo corr)
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/prototype/' cpa_id '/QA_peak_topo/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Logic for spatial combo
both_ns_idx   = ~any([peak_idx, topo_corr_idx],2);
pk_only_idx   = peak_idx & ~topo_corr_idx;
topo_only_idx = ~peak_idx & topo_corr_idx;
pk_mrkr   = 'd';
topo_mrkr = 's';
ns_mrkr   = 'o';

% Combined spatial (peak elec vs. topo corr)
fig_name = [SBJ '_proto_' cpa_id '_QA_peak_topo'];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
hold on;
scatter(n_elec_match(both_ns_idx),topo_corrs(both_ns_idx),ns_sz,'k',ns_mrkr);
scatter(n_elec_match(pk_only_idx),topo_corrs(pk_only_idx),sig_sz,'r',pk_mrkr);
scatter(n_elec_match(topo_only_idx),topo_corrs(topo_only_idx),sig_sz,'r',topo_mrkr);
scatter(n_elec_match(final_ics),topo_corrs(final_ics),sig_sz,'r','filled');
xlim([0 cpa.n_max_elec]);
ylim([-1 1]);
line([cpa.min_elec_match cpa.min_elec_match], ylim, 'Color','r');
% line([cpa.ic_rank_max cpa.ic_rank_max],ylim,'Color','r');
xlabel('# elec match');
ylabel('topo corr');
title([num2str(sum(space_ic_idx)) '/' num2str(numel(final_ics)) '; peak: ' ...
    num2str(sum(peak_idx)) '; corr: ' num2str(sum(topo_corr_idx))]);
set(gca,'FontSize',14);

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Combined temporal-spatial (Peak vs. temporal)
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/prototype/' cpa_id '/QA_peak_time/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

fig_name = [SBJ '_proto_' cpa_id '_QA_peak_time'];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);

% Logic for spatial combo
both_ns_idx   = ~any([peak_idx, time_ic_idx],2);
pk_only_idx   = peak_idx & ~time_ic_idx;
time_only_idx = ~peak_idx & time_ic_idx;
time_mrkr = '*';

% Plot scatter of both metrics
ylim([0 diff(cpa.time_win)]);
hold on;
scatter(n_elec_match(both_ns_idx),sig_len(both_ns_idx),ns_sz,'k',ns_mrkr);
scatter(n_elec_match(pk_only_idx),sig_len(pk_only_idx),sig_sz,'r',pk_mrkr);
scatter(n_elec_match(time_only_idx),sig_len(time_only_idx),sig_sz,'r',time_mrkr);
scatter(n_elec_match(final_ics),sig_len(final_ics),sig_sz,'r','filled');
xlim([0 cpa.n_max_elec]);
ylim([0 diff(cpa.time_win)]);
line([cpa.min_elec_match cpa.min_elec_match], ylim, 'Color','r');
line(xlim, [cpa.min_sig_len cpa.min_sig_len],'Color','r');
xlabel('# elec match');
ylabel('sig len (s)');
title(['peak: ' num2str(sum(peak_idx)) '; time: ' num2str(sum(time_ic_idx))]);
set(gca,'FontSize',14);

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Save Data
clean_data_fname = [SBJ_vars.dirs.proc SBJ '_' cpa_id '_' proc_id '_prototype.mat'];
fprintf('Saving %s\n',clean_data_fname);
save(clean_data_fname, '-v7.3', 'final_ics','clean_ica','sig_wins','sig_len',...
    'sig_len_chunk','top_elecs','topo_corrs','topo_pvals');

end
