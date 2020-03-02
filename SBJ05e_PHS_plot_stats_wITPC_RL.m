function SBJ05e_PHS_plot_stats_wITPC_RL(SBJs,proc_id,an_id,stat_id,save_fig,varargin)
%% Plot group TFRs per regressor: z-scored weighted ITPC outlined by significance
%   Only for single channel right now...
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
        elseif strcmp(varargin{v},'plt_id') && ischar(varargin{v+1})
            plt_id = varargin{v+1};
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
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~an.itpc; error('why run this without ITPC an_vars?'); end

stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.an_style,'wITPC'); error('stat_id not using wITPC!'); end
% plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
% eval(plt_vars_cmd);

% Select Conditions of Interest
[reg_lab, reg_names, reg_colors, reg_styles]  = fn_regressor_label_styles(st.model_lab);

%% Load Stats
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' stat_id '_' an_id '.mat'],'SBJs');
if ~all(strcmp(SBJs,tmp.SBJs))
    fprintf(2,'Loaded SBJs: %s\n',strjoin(tmp.SBJs,', '));
    error('Not all SBJs match input SBJ list!');
end
warning('WARNING: Assuming same prdm_vars for all SBJ to get event timing!');
prdm_vars = load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/03_events/' SBJs{1} '_prdm_vars.mat']);

load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' stat_id '_' an_id '.mat'],'zwitpc','qvals');

%% Load TFR for axes
load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/04_proc/' SBJs{1} '_' proc_id '_' an_id '.mat']);
if numel(tfr.label) > 1; error('only ready for one channel right now!'); end

% Select time and trials of interest
cfgs = []; cfgs.latency = st.stat_lim;
st_tfr = ft_selectdata(cfgs, tfr);
st_time_vec = st_tfr.time;
fois = st_tfr.freq;

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' an_id '/' stat_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(st_tfr.label)
    % Get color lims per condition
    data_mat = squeeze(mean(zwitpc,1));
    clim = [min(data_mat(:)) max(data_mat(:))];
    
    % Get significance mask
    ns_alpha = 0.4;
    sig_mask = ones(size(qvals))*ns_alpha;
    sig_mask(qvals<=st.alpha) = 1;

%     tick_ix = 1:3:numel(fois);
%     yticklab = cell(size(tick_ix));
%     for f = 1:numel(tick_ix)
%         yticklab{f} = num2str(st_tfr.freq(tick_ix(f)),'%.1f');
%     end
    
    %% Create plot
    fig_name = ['GRP_' stat_id '_' an_id '_' st_tfr.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Beta Plots
    %cfgplt = []; cfgplt.zlim = clim;
    axes = gobjects([numel(reg_lab) 1]);
    [num_rc,~] = fn_num_subplots(numel(reg_lab));
    for reg_ix = 1:numel(reg_lab)
        subplot(num_rc(1),num_rc(2),reg_ix);
        axes(reg_ix) = gca; hold on;
        
        im = imagesc(st_time_vec, fois, squeeze(data_mat(reg_ix,:,:)),clim);
        im.AlphaData = squeeze(sig_mask(reg_ix,:,:));
        set(axes(reg_ix),'YDir','normal');
%         set(axes(reg_ix),'YTick',tick_ix);
%         set(axes(reg_ix),'YTickLabels',yticklab);
        set(axes(reg_ix),'YLim',[min(fois) max(fois)]);
        set(axes(reg_ix),'XLim',[min(st_time_vec) max(st_time_vec)]);
        set(axes(reg_ix),'XTick',[min(st_time_vec):0.1:max(st_time_vec)]);
        %ft_singleplotTFR(cfgplt, tfr_avg{cond_ix});
        title([st_tfr.label{ch_ix} ': z-wITPC ' reg_lab{reg_ix}]);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar('northoutside');
        set(axes(reg_ix),'FontSize',16);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.png'];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end

end
