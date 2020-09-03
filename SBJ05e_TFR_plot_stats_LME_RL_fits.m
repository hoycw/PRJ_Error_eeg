function SBJ05e_TFR_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_id,stat_id,save_fig,varargin)
%% Plot group model coefficient TFR matrix per regressor + R2 plot
%   Betas are outlined by significance
%   Maximum effect size per regressor is marked with scatter point
%   Only for single channel
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
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
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select Conditions of Interest
[reg_lab, reg_names, ~, ~]  = fn_regressor_label_styles(st.model_lab);

%% Load Stats
% Check stats are run on correct SBJ group
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'SBJs');
if ~all(strcmp(SBJs,tmp.SBJs))
    fprintf(2,'Loaded SBJs: %s\n',strjoin(tmp.SBJs,', '));
    error('Not all SBJs match input SBJ list!');
end

% Load TFR stats
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'lme','qvals');

%% Load Example TFR for Axes
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
    % Get color limits per regressor
    beta_mat = zeros([numel(reg_lab) numel(fois) numel(st_time_vec)]);
    r2_mat   = zeros([numel(fois) numel(st_time_vec)]);
    for f_ix = 1:numel(fois)
        for t_ix = 1:numel(st_time_vec)
            beta_mat(:,f_ix,t_ix) = lme{f_ix,t_ix}.Coefficients.Estimate(2:end);
            r2_mat(f_ix,t_ix)     = lme{f_ix,t_ix}.Rsquared.Adjusted;
        end
    end
    
    % Find color limits across regressors
    clim = [-max(abs(beta_mat(:))) max(abs(beta_mat(:)))];
    % Fix color limits across analyses (channels) for final paper plots
    if strcmp(SBJ_id,'goodall')
        clim = [-1.5 1.5]; % Max beta is -1.4 for sRPE in Fz
    end
    
    % Get significance mask
    ns_alpha = 0.4;
    sig_mask = ones([numel(reg_lab) numel(fois) numel(st_time_vec)])*ns_alpha;
    sig_mask(qvals<=st.alpha) = 1;
    
    % Find max beta points for markers
    beta_pks = zeros([numel(reg_lab) 2]);
    beta_pk_f_ix = zeros([numel(reg_lab) 2]);
    beta_pk_t_ix = zeros([numel(reg_lab) 2]);
    minmax_ix = zeros(size(reg_lab));
    for reg_ix = 1:numel(reg_lab)
        for f = 1:numel(fois)
            % Find minimum beta
            if min(beta_mat(reg_ix,f,:)) < beta_pks(reg_ix,1)
                [~, beta_pk_t_ix(reg_ix,1)] = min(beta_mat(reg_ix,f,:));
                beta_pk_f_ix(reg_ix,1) = f;
                beta_pks(reg_ix,1) = beta_mat(reg_ix,beta_pk_f_ix(reg_ix,1),beta_pk_t_ix(reg_ix,1));
            end
            % Find maximum beta
            if max(beta_mat(reg_ix,f,:)) > beta_pks(reg_ix,2)
                [~, beta_pk_t_ix(reg_ix,2)] = max(beta_mat(reg_ix,f,:));
                beta_pk_f_ix(reg_ix,2) = f;
                beta_pks(reg_ix,2) = beta_mat(reg_ix,beta_pk_f_ix(reg_ix,2),beta_pk_t_ix(reg_ix,2));
            end
        end
        
        % Take beta with largest absolute value
        if abs(beta_pks(reg_ix,1)) > abs(beta_pks(reg_ix,2))
            minmax_ix(reg_ix) = 1;
        else
            minmax_ix(reg_ix) = 2;
        end
    end
    
    %% Create plot
    fig_name = [SBJ_id '_' stat_id '_' an_id '_' st_tfr.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Beta Plots
    axes = gobjects([numel(reg_lab)+1 1]);
    [num_rc,~] = fn_num_subplots(numel(reg_lab)+1);
    for reg_ix = 1:numel(reg_lab)
        subplot(num_rc(1),num_rc(2),reg_ix);
        axes(reg_ix) = gca; hold on;
        
        % Plot Beta Matrix
        im = imagesc(st_time_vec, fois, squeeze(beta_mat(reg_ix,:,:)),clim);
        im.AlphaData = squeeze(sig_mask(reg_ix,:,:));
        
        % Plot max/min point
        scatter(st_time_vec(beta_pk_t_ix(reg_ix,minmax_ix(reg_ix))),...
            fois(beta_pk_f_ix(reg_ix,minmax_ix(reg_ix))), 100,'r','*');
        
        % Plot Properties
        set(axes(reg_ix),'YDir','normal');
        set(axes(reg_ix),'YLim',[min(fois) max(fois)]);
        set(axes(reg_ix),'XLim',[min(st_time_vec) max(st_time_vec)]);
        set(axes(reg_ix),'XTick',[min(st_time_vec):0.1:max(st_time_vec)]);
        title([st_tfr.label{ch_ix} ': ' reg_names{reg_ix} ' Beta']);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar('northoutside');
        set(axes(reg_ix),'FontSize',16);
    end
    
    % R2 Plot
    subplot(num_rc(1),num_rc(2),numel(reg_lab)+1);
    axes(numel(reg_lab)+1) = gca; hold on;
    
    % Plot R2 Matrix
    imagesc(st_time_vec, fois, r2_mat);
    
    % Axes and labels
    set(axes(numel(reg_lab)+1),'YDir','normal');
    set(axes(numel(reg_lab)+1),'YLim',[min(fois) max(fois)]);
    set(axes(numel(reg_lab)+1),'XLim',[min(st_time_vec) max(st_time_vec)]);
    set(axes(numel(reg_lab)+1),'XTick',[min(st_time_vec):0.1:max(st_time_vec)]);
    title([st_tfr.label{ch_ix} ': R2 (n=' num2str(numel(SBJs)) ')']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar('northoutside');
    set(gca,'FontSize',16);
    
    % Report min and max beta points and stats
    for reg_ix = 1:numel(reg_lab)
        fprintf('min %s = %.06f at %.03f s, %.02f Hz; p = %.10f\n',reg_lab{reg_ix},...
            beta_pks(reg_ix,1),st_time_vec(beta_pk_t_ix(reg_ix,1)),...
            fois(beta_pk_f_ix(reg_ix,1)),qvals(reg_ix,beta_pk_f_ix(reg_ix,1),beta_pk_t_ix(reg_ix,1)));
        fprintf('max %s = %.06f at %.03f s, %.02f Hz; p = %.10f\n',reg_lab{reg_ix},...
            beta_pks(reg_ix,2),st_time_vec(beta_pk_t_ix(reg_ix,2)),...
            fois(beta_pk_f_ix(reg_ix,2)),qvals(reg_ix,beta_pk_f_ix(reg_ix,2),beta_pk_t_ix(reg_ix,2)));
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end

end
