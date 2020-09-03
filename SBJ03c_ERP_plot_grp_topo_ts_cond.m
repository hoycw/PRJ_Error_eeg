function SBJ03c_ERP_plot_grp_topo_ts_cond(SBJ_id,conditions,proc_id,an_id,stat_ids,save_fig,varargin)
%% Plot group ERP topographies per condition across multiple windows (dynamics)
%   Plotting evoked activity in time window identified as important via
%   stats, e.g., based on times used in SBJ04c_ERP_grp_stats_LME_RL
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   conditions [str] - group of condition labels to segregate trials
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_ids [cell array] - IDs of statistical analyses that provide peak times
%   save_fig [0/1] - binary flag to save figure
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   saves figure

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
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

%% Load Results
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
if ~strcmp(an.ROI{1},'all'); error('run this only for topos with all channels!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select conditions (and trials)
[cond_lab, cond_names, ~, ~, ~] = fn_condition_label_styles(conditions);

%% Load Plotting Windows from stat_id list
sts      = cell(size(stat_ids));
pk_times = zeros(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_ids{st_ix} '_vars.m'];
    eval(stat_vars_cmd);
    sts{st_ix} = st;
    if ~strcmp(st.measure,'mean'); error('this script is for mean window only!');end
    
    % Load Peak Time
    load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_ids{st_ix} '_' an_id '.mat'],'ch_list','reg_pk_time');
    pk_times(st_ix) = reg_pk_time;
    
    % Ensure all stat params are equal except peak center
    if st_ix>1
        fnames = fieldnames(sts{st_ix});
        for f_ix = 1:numel(fnames)
            if ~any(strcmp(fnames{f_ix},{'pk_center','pk_reg_id','pk_an_id'}))
                if ischar(sts{st_ix}.(fnames{f_ix}))
                    if ~strcmp(sts{1}.(fnames{f_ix}), sts{st_ix}.(fnames{f_ix}))
                        error(['st.' fnames{f_ix} ' not the same!']);
                    end
                elseif isnumeric(sts{st_ix}.(fnames{f_ix}))
                    if any(sts{1}.(fnames{f_ix}) ~= sts{st_ix}.(fnames{f_ix}))
                        error(['st.' fnames{f_ix} ' not the same!']);
                    end
                elseif iscell(sts{st_ix}.(fnames{f_ix}))
                    if numel(sts{st_ix}.(fnames{f_ix}))>1; error('not ready for cells n > 1'); end
                    if any(sts{1}.(fnames{f_ix}){1} ~= sts{st_ix}.(fnames{f_ix}){1})
                        error(['st.' fnames{f_ix} ' not the same!']);
                    end
                else
                    error(['Unknown class(' fnames{f_ix} ')']);
                end
            end
        end
    end
    clear st stat_vars_cmd
end

%% Load data
er_avg   = cell([numel(cond_lab) numel(stat_ids) numel(SBJs)]);
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',an_id,'.mat'],'roi');
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/03_events/',SBJs{s},'_behav_',proc_id,'_final.mat'],'bhv');
    
    % Separate out each condition
    cfg_er = [];
    cfg_avgtime = []; cfg_avgtime.avgovertime = 'yes';
    for cond_ix = 1:numel(cond_lab)
        cond_idx = fn_condition_index(cond_lab(cond_ix),bhv);
        cfg_er.trials = find(cond_idx);
        for st_ix = 1:numel(stat_ids)
            % Average ERP within stat_id window
            cfg_er.latency = sts{st_ix}.stat_lim + pk_times(st_ix);
            er_avg{cond_ix,st_ix,s} = ft_timelockanalysis(cfg_er,roi);
            er_avg{cond_ix,st_ix,s} = ft_selectdata(cfg_avgtime,er_avg{cond_ix,st_ix,s});
        end
    end
    clear roi bhv
end

%% Compute group Grand Average ERPs for plotting
er_grp = cell([numel(cond_lab) numel(stat_ids)]);
clim   = [nan nan];
for cond_ix = 1:numel(cond_lab)
    for st_ix = 1:numel(stat_ids)
        er_grp{cond_ix,st_ix} = ft_timelockgrandaverage([], er_avg{cond_ix,st_ix,:});
        
        % Get color limits for plotting
        clim(1) = min([clim(1); er_avg{cond_ix,st_ix}.avg]);
        clim(2) = max([clim(2); er_avg{cond_ix,st_ix}.avg]);
    end
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' strjoin(stat_ids,'-') '/' conditions '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create plot
fig_name = [SBJ_id '_topo_ts_' conditions '_' an_id];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);

%% Create a figure for each channel
% Set up default plotting params
axes = gobjects([numel(cond_lab) numel(stat_ids)]);
cfgp = [];
cfgp.zlim      = clim;
cfgp.layout    = 'biosemi64.lay';
cfgp.comment   = 'no';
cfgp.highlight = 'off';
cfgp.parameter = 'avg';
cfgp.colormap  = 'jet';

% Plot topo for each time and condition
plot_ix = 0;
for cond_ix = 1:numel(cond_lab)
    for st_ix = 1:numel(stat_ids)
        plot_ix = plot_ix + 1;
        subplot(numel(cond_lab),numel(stat_ids),plot_ix);
        axes(cond_ix,st_ix) = gca; hold on;
        
        % Only plot color bar on last column
        if st_ix==numel(stat_ids)
            cfgp.colorbar = 'yes';
        else
            cfgp.colorbar = 'no';
        end
        
        % Plot ERP topo
        cfgp.xlim   = [pk_times(st_ix) pk_times(st_ix)];
        tmp = ft_topoplotER(cfgp, er_grp{cond_ix,st_ix});
        title([cond_names{cond_ix} ': ' num2str(pk_times(st_ix),'%.3f') ' s']);
        axis tight
    end
end

%% Save figure
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
