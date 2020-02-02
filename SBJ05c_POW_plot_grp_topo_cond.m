function SBJ05c_POW_plot_grp_topo_cond(SBJs,conditions,proc_id,an_id,plt_id,save_fig,varargin)
%% Plot single-band power topography per condition for single window across group
% INPUTS:
%   conditions [str] - group of condition labels to segregate trials

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
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
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select conditions (and trials)
[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);

%% Load data
er_avg = cell([numel(cond_lab) numel(SBJs)]);
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/03_events/',SBJs{s},'_behav_',proc_id,'_final.mat'],'bhv');
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',proc_id,'_',an_id,'.mat'],'tfr');
    if numel(tfr.freq)>1; error('TFR is not averaged over frequencies!'); end
    if numel(tfr.label)==1; error([an_id ' has only one channel!']); end
    
    % Separate out each condition
    cfg_er = [];
    for cond_ix = 1:numel(cond_lab)
        cond_idx = fn_condition_index(cond_lab(cond_ix),bhv);
        cfg_er.trials = find(cond_idx);
        er_avg{cond_ix,s} = ft_timelockanalysis(cfg_er,tfr);
        
        % Remove frequency from label so topo plot doesn't error
        for e_ix = 1:numel(tfr.label)
            er_avg{cond_ix,s}.label{e_ix} = ...
                er_avg{cond_ix,s}.label{e_ix}(1:strfind(er_avg{cond_ix,s}.label{e_ix},'@')-1);
        end
    end
    clear tfr bhv
end

%% Average ERPs for plotting
cfg_gavg = [];
er_grp = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    er_grp{cond_ix} = ft_timelockgrandaverage(cfg_gavg, er_avg{cond_ix,:});
end

% Get color limits
cfgat = [];
cfgat.latency = plt.plt_lim;
cfgat.avgovertime = 'yes';
clim = [0 0];
for cond_ix = 1:numel(cond_lab)
    tmp = ft_selectdata(cfgat,er_grp{cond_ix});
    clim = [min([clim(1) min(tmp.avg)]) max([clim(2) max(tmp.avg)])];
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' an_id '/' conditions '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create plot
fig_name = ['GRP_' conditions '_' an_id '_' plt_id];
[num_rc,~] = fn_num_subplots(numel(cond_lab));
if num_rc(1)>1; fig_height = 0.5; else, fig_height = 0.3; end
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 fig_height],'Visible',fig_vis);   %this size is for single plots

%% Create a figure for each channel
axes = gobjects(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    subplot(num_rc(1),num_rc(2),cond_ix);
    axes(cond_ix) = gca; hold on;
    
    cfgp = [];
    cfgp.xlim     = plt.plt_lim;
    cfgp.zlim     = clim;
    cfgp.layout   = 'biosemi64.lay';
    cfgp.colorbar = 'yes';
    cfgp.comment  = 'no';
    tmp = ft_topoplotER(cfgp, er_grp{cond_ix});
    title(cond_lab{cond_ix});
    tmp = caxis;
    clim = [min([clim(1) tmp(1)]) max([clim(2) tmp(2)])]; 
    axis tight
end

%% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
