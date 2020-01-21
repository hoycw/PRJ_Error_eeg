function SBJ06a_CPA_prototype_selection(SBJ, proc_id, cpa_id, plt_id,save_fig,varargin)
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
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load the data
%loaded data from after SBJ02a --> already cleaned and trial segmented
load([SBJ_vars.dirs.preproc SBJ '_preproc_eeg_full_ft.mat']);
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat']); %chose 02a - ica before rejection!
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles('Odd'); % maybe change this so not hardcoded
cond_idx = fn_condition_index(cond_lab, bhv);

% Create contrast: (Unexpected - Expected) for each outcome
[diff_lab, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles(cpa.diff_id);
if numel(diff_lab)>1; error('Too many condition contrasts!'); end

%% Compute plotting data
trials = cell(size(cond_lab));
means  = NaN([numel(cond_lab) numel(ica.label) numel(ica.time{1})]);
sems   = NaN([numel(cond_lab) numel(ica.label) numel(ica.time{1})]);
for cond_ix = 1:numel(cond_lab)
    % Select trials for plotting
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(ica.label) numel(cond_trial_ix) numel(ica.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,t_ix,:) = ica.trial{cond_trial_ix(t_ix)};
    end
    
    % Compute mean and variance
    means(cond_ix,:,:) = squeeze(mean(trials{cond_ix},2));
    sems(cond_ix,:,:) = squeeze(std(trials{cond_ix},[],2))./sqrt(size(trials{cond_ix},2))';
end

%% Temporal Criteria: Condition Stats
% Select stats epoch
cfg = [];
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

final_ics = find(all([time_ic_idx, space_ic_idx, var_ic_idx],2));
fprintf('%s: %d / %d components selected!\n',SBJ,numel(final_ics),numel(ica.label));
if sum(final_ics)<1
    error('No ICs found that match all criteria!');
end

%% PLOT AND SAVE
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/prototype/' cpa_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

for f_ix = 1:numel(final_ics)
    comp_ix = final_ics(f_ix);
    fig_name = [SBJ '_' proc_id '_' cpa_id '_IC' num2str(comp_ix)];
    figure('Name', fig_name, 'units','normalized',...
        'outerposition',[0 0 0.6 0.5], 'Visible', fig_vis);
    axes = subplot(1,3,[1 2]); hold on;
    
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+1 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(ica.time{1}, means(cond_ix, comp_ix, :), sems(cond_ix, comp_ix, :),...
            {'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},plt.errbar_alpha);
        hold on
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    % Plot Significance
    [sig_lims] = fn_find_chunks(sig_wins(comp_ix,:));
    sig_lims(squeeze(sig_wins(comp_ix,sig_lims(:,1)))==0,:) = [];
    
    data_lim = [min(min(means(:,comp_ix,:)-sems(:,comp_ix,:))) max(max(means(:,comp_ix,:)+sems(:,comp_ix,:)))];
    for sig_ix = 1:size(sig_lims,1)
        sig_times = st_ica.time{1}(sig_lims(sig_ix,1):sig_lims(sig_ix,2));
        sig_y = data_lim(1) + data_lim(1)*plt.sig_loc_factor;
        sig_line = line(sig_times,repmat(sig_y,size(sig_times)),...
            'LineWidth',plt.sig_width,'Color',diff_colors{1});
        if sig_ix==1
            main_lines(end) = sig_line;
        end
    end
    
    line([0 0],ylim,'Color','k');
    % Axes and Labels
    axes(1).YLabel.String = 'uV';
    axes(1).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(1).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(1).XLabel.String = 'Time (s)';
    title(['IC#' num2str(comp_ix) '- sig len (s): '...
        num2str(sig_len(comp_ix)/ica.fsample) ' (' num2str(100*sig_len(comp_ix)/numel(st_ica.time{1}),'%.1f') '%)']);
    leg_lab = [cond_lab diff_lab(1)];%'F' 
    if plt.legend
       legend(main_lines,leg_lab,'Location',plt.legend_loc);
    end
    set(gca,'FontSize',16);
    
    % Plot Topo
    subplot(1,3,3); hold on;
    cfgp = [];
    cfgp.component = comp_ix;
%     cfgp.highlightchannel = top_elecs(comp_ix,:);
%     cfgp.highlightsymbol = '*';
    cfgp.layout    = 'biosemi64.lay';
    cfgp.comment   = 'no';
    ft_topoplotIC(cfgp, ica);
    title(['Peak Elecs: ' strjoin(top_elecs(comp_ix,:),',')]); 
    set(gca,'FontSize',16);

    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Save Data
clean_data_fname = [SBJ_vars.dirs.proc SBJ '_' cpa_id '_' proc_id '_prototype.mat'];
save(clean_data_fname, '-v7.3', 'final_ics');

end
