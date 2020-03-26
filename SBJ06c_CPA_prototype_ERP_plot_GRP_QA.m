function SBJ06c_CPA_prototype_ERP_plot_GRP_QA(SBJ_id,conditions,proc_id,cpa_id,save_fig,varargin)
error('added this to SBJ06a for single SBJ and all ICs, never wrote this group level script');
%% Plot QA for prototype selection across group
%   Plot significant time vs. IC #
% INPUTS:
%   conditions [str] - group of condition labels to segregate trials

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
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
if ~strcmp(conditions,'Odd'); error('Why run prototype ERPs if not Odd?'); end

%% Load Results
cpa_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' cpa_id '_vars.m'];
eval(cpa_vars_cmd);

% Select SBJs
sbj_file = fopen([root_dir 'PRJ_Error_EEG/scripts/SBJ_lists/' SBJ_id '.sbj']);
tmp = textscan(sbj_file,'%s');
fclose(sbj_file);
SBJs = tmp{1}; clear tmp;

[cond_lab, cond_names, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);

%% Load example data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m'];
eval(SBJ_vars_cmd);
tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/04_proc/',...
    SBJs{1} '_' cpa_id '_' proc_id '_prototype.mat'],'clean_ica');
cfg = []; cfg.latency = cpa.time_win;
st_ica = ft_selectdata(cfg,tmp.clean_ica);
clear tmp

%% Load all data
ic_list = cell(size(SBJs));
sig_len = cell(size(SBJs));
space_metric = cell(size(SBJs));
for s = 1:length(SBJs)
    % Load SBJ data
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    load([SBJ_vars.dirs.proc SBJs{s} '_' cpa_id '_' proc_id '_prototype.mat']);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    ic_list{s} = final_ics;
    
    % Compile Temporal Criteria
    sig_len{s} = nan(size(final_ics));
    for f_ix = 1:numel(final_ics)
        comp_ix = final_ics(f_ix);
        [sig_lims] = fn_find_chunks(sig_wins(comp_ix,:));
        sig_lims(squeeze(sig_wins(comp_ix,sig_lims(:,1)))==0,:) = [];
        for seg_ix = 1:size(sig_lims,1)
            sig_len{s}(f_ix) = sum(diff(sig_lims,1,2)+1)./clean_ica.fsample;
        end
    end
    
    % Compile Spatial Criteria
    space_metric{s} = nan(size(final_ics));
    if strcmp(cpa.elec_method,'peak')
        for f_ix = 1:numel(final_ics)
            comp_ix = final_ics(f_ix);
            space_metric{s}(f_ix) = numel(intersect(top_elecs(comp_ix,:),cpa.elec_list));
        end
    elseif strcmp(cpa.elec_method,'topo_corr')
        for f_ix = 1:numel(final_ics)
            comp_ix = final_ics(f_ix);
            space_metric{s}(f_ix) = topo_corrs(comp_ix);
        end
    end
     
    clear tmp SBJ_vars bhv clean_ica final_ics sig_wins top_elecs topo_corrs
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/CPA/prototype/' cpa_id '/' conditions '/' an_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Create plot
    fig_name = [SBJ_id '_proto_' cpa_id '_' conditions '_' an_id '_' ch_list{ch_ix}];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);   %this size is for single plots
%     [plot_rc,~] = fn_num_subplots(numel(w2.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     subplot(plot_rc(1),plot_rc(2),ch_ix);
    ax = gca; hold on;
    
    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(an.event_type) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(time_vec, plt_means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_lines(numel(cond_lab)+evnt_ix) = line(...
            [evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end
    leg_lab = [cond_lab plt.evnt_lab];
    
    % Axes and Labels
    ax.YLabel.String = 'uV';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    ax.Title.String  = [ch_list{ch_ix} ' (n=' num2str(numel(SBJs)) ')'];
    set(ax,'FontSize',16');
    if plt.legend
        legend(main_lines,leg_lab{:},'Location',plt.legend_loc);
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        % Ensure vector graphics if saving
        if any(strcmp(fig_ftype,{'svg','eps'}))
            set(gcf, 'Renderer', 'painters');
        end
        saveas(gcf,fig_fname);
    end
end

end
