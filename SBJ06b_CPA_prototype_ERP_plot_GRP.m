function SBJ06b_CPA_prototype_ERP_plot_GRP(SBJ_id, proc_id, cpa_id, an_id, plt_id, save_fig,varargin)
error('use broken out reconstruct, ERP, plot scripts');
%% Candidate-Prototype Analysis (CPA): Prototype ERP Plot
%   Reconstruct oddball data based on final prototype ICs
%   Compute ERPs, and plot at group level

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
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
cpa_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' cpa_id '_vars.m'];
eval(cpa_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
sbj_file = fopen([root_dir 'PRJ_Error_EEG/scripts/SBJ_lists/' SBJ_id '.sbj']);
tmp = textscan(sbj_file,'%s');
fclose(sbj_file);
SBJs = tmp{1}; clear tmp;

[cond_lab, cond_names, cond_colors, cond_styles, ~] = fn_condition_label_styles('Odd');

% Create contrast: (Unexpected - Expected) for each outcome
[diff_lab, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles(cpa.diff_id);
if numel(diff_lab)>1; error('Too many condition contrasts!'); end
    
%% Load the data
% Load example data
load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/04_proc/' SBJs{1} '_' ...
    cpa_id '_' proc_id '_prototype.mat'],'clean_ica');
plot_time = clean_ica.time{1};
% Select stats epoch
cfg = []; cfg.latency = cpa.time_win;
st_ica = ft_selectdata(cfg,clean_ica);

% Load EEG20 (full cap, no bad chan) for topo channel order
% tmp = load([root_dir 'PRJ_Error_eeg/data/EEG01/02_preproc/EEG01_eeg_full_ft_final.mat'],'clean_trials');
% ch_list = tmp.clean_trials.label; clear tmp

% topo   = NaN([numel(SBJs) numel(ch_list)]);
means  = NaN([numel(cond_lab) numel(SBJs) numel(clean_ica.time{1})]);
sems   = NaN([numel(cond_lab) numel(SBJs) numel(clean_ica.time{1})]);
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    proto_fname = [SBJ_vars.dirs.proc SBJs{s} '_' cpa_id '_' proc_id '_prototype.mat'];
    if exist(proto_fname,'file')
        load(proto_fname);
        load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
        cond_idx = fn_condition_index(cond_lab, bhv);
        
        % Reconstruction
        cfg = [];
        cfg.component = setdiff(1:numel(clean_ica.label),final_ics);
        recon = ft_rejectcomponent(cfg, clean_ica);
        
        % Repair Bad Channels
        cfg = [];
        cfg.method         = 'average';
        cfg.missingchannel = SBJ_vars.ch_lab.bad(:); % not in data (excluded from ica)
        cfg.layout         = 'biosemi64.lay';
        
        cfgn = [];
        cfgn.channel = 'all';
        cfgn.layout  = 'biosemi64.lay';
        cfgn.method  = 'template';
        cfg.neighbours = ft_prepare_neighbours(cfgn);
        full_recon = ft_channelrepair(cfg, recon);
        
        Average topo
        topo_sbj = nan([size(topo,2) 1]);
        for ch_ix = 1:numel(clean_ica.topolabel)
            full_ch_ix = find(strcmp(ch_list,clean_ica.topolabel{ch_ix}));
            topo_sbj(full_ch_ix) = mean(clean_ica.topo(ch_ix,final_ics));
        end
        topo(s,:) = sbj_topo;
        
        % Compute plotting data
        for cond_ix = 1:numel(cond_lab)
            % Select trials for plotting
            cond_trial_ix = find(cond_idx==cond_ix);
            cond_trials = nan([numel(final_ics) numel(cond_trial_ix) numel(clean_ica.time{1})]);
            for t_ix = 1:numel(cond_trial_ix)
                cond_trials(:,t_ix,:) = clean_ica.trial{cond_trial_ix(t_ix)}(final_ics,:);
            end
            
            % Compute mean and variance
            means(cond_ix,s,:) = squeeze(mean(cond_trials,2));
            sems(cond_ix,s,:) = squeeze(std(cond_trials,[],2))./sqrt(size(cond_trials,2))';
        end
    else
        fprintf(2,'\tWarning: %s does not exist, skipping %s\n',proto_fname,SBJs{s});
    end
    
    clear SBJ_vars full_recon final_ics clean_ica cond_idx cond_trials cond_trial_ix
end

%% Combien across SBJ
plot_mean = mean(means,2);
plot_sem  = mean(sems,2);

%% Plot data
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/prototype/' cpa_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

for f_ix = 1:numel(final_ics)
    comp_ix = final_ics(f_ix);
    fig_name = [SBJ_id '_' proc_id '_' cpa_id];
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
