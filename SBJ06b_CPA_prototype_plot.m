function SBJ06b_CPA_prototype_plot(SBJ, proc_id, cpa_id, plt_id,save_fig,varargin)
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
load([SBJ_vars.dirs.proc SBJ '_' cpa_id '_' proc_id '_prototype.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

[cond_lab, ~, cond_colors, cond_styles, ~] = fn_condition_label_styles('Odd'); % maybe change this so not hardcoded
cond_idx = fn_condition_index(cond_lab, bhv);

% Create contrast: (Unexpected - Expected) for each outcome
[diff_lab, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles(cpa.diff_id);
if numel(diff_lab)>1; error('Too many condition contrasts!'); end

%% Compute plotting data
trials = cell(size(cond_lab));
means  = NaN([numel(cond_lab) numel(clean_ica.label) numel(clean_ica.time{1})]);
sems   = NaN([numel(cond_lab) numel(clean_ica.label) numel(clean_ica.time{1})]);
for cond_ix = 1:numel(cond_lab)
    % Select trials for plotting
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(clean_ica.label) numel(cond_trial_ix) numel(clean_ica.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,t_ix,:) = clean_ica.trial{cond_trial_ix(t_ix)};
    end
    
    % Compute mean and variance
    means(cond_ix,:,:) = squeeze(mean(trials{cond_ix},2));
    sems(cond_ix,:,:) = squeeze(std(trials{cond_ix},[],2))./sqrt(size(trials{cond_ix},2))';
end

% Select stats epoch
cfg = [];
cfg.latency = cpa.time_win;
st_ica = ft_selectdata(cfg,clean_ica);

%% Plot data
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/prototype/' cpa_id '/final_ics/'];
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
        ebars{cond_ix} = shadedErrorBar(clean_ica.time{1}, means(cond_ix, comp_ix, :), sems(cond_ix, comp_ix, :),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        hold on
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    % Compute summary metrics (% sig, min_sig_len)
    %sig_perc = sum(sig_wins(comp_ix,:)) / size(sig_wins,2);
    [sig_lims] = fn_find_chunks(sig_wins(comp_ix,:));
    sig_lims(squeeze(sig_wins(comp_ix,sig_lims(:,1)))==0,:) = [];
    if ~isempty(sig_lims)
        sig_len = sum(diff(sig_lims,1,2)+1);
    end
    
    % Plot Significance
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
        num2str(sig_len/clean_ica.fsample) ' ('...
        num2str(100*sig_len/numel(st_ica.time{1}),'%.1f') '%)']);
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
    ft_topoplotIC(cfgp, clean_ica);
    if strcmp(cpa.elec_method,'peak')
        space_str = strjoin(top_elecs(comp_ix,:),',');
    elseif strcmp(cpa.elec_method,'topo_corr')
        space_str = ['r=' num2str(topo_corrs(comp_ix),'%.3f')];
    end
    title(['Peak Elecs: ' space_str]); 
    set(gca,'FontSize',16);

    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
