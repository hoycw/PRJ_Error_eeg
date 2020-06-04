function SBJ05e_TFR_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_id,stat_id,save_fig,varargin)
%% Plot group TFRs per regressor: betas outlined by significance; + R2 plot
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
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
% plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
% eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select Conditions of Interest
[reg_lab, ~, reg_colors, reg_styles]  = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, cond_colors, cond_styles, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Stats
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'SBJs');
if ~all(strcmp(SBJs,tmp.SBJs))
    fprintf(2,'Loaded SBJs: %s\n',strjoin(tmp.SBJs,', '));
    error('Not all SBJs match input SBJ list!');
end
warning('WARNING: Assuming same prdm_vars for all SBJ to get event timing!');
prdm_vars = load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/03_events/' SBJs{1} '_prdm_vars.mat']);

load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'lme','qvals');

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
    beta_mat = zeros([numel(reg_lab) numel(fois) numel(st_time_vec)]);
    r2_mat   = zeros([numel(fois) numel(st_time_vec)]);
    for f_ix = 1:numel(fois)
        for t_ix = 1:numel(st_time_vec)
            beta_mat(:,f_ix,t_ix) = lme{f_ix,t_ix}.Coefficients.Estimate(2:end);
            r2_mat(f_ix,t_ix)     = lme{f_ix,t_ix}.Rsquared.Adjusted;
        end
    end
    clim = [-max(abs(beta_mat(:))) max(abs(beta_mat(:)))];
    
    % Get significance mask
    ns_alpha = 0.4;
    sig_mask = ones([numel(reg_lab) numel(fois) numel(st_time_vec)])*ns_alpha;
    sig_mask(qvals<=st.alpha) = 1;
    
%     tick_ix = 1:3:numel(fois);
%     yticklab = cell(size(tick_ix));
%     for f = 1:numel(tick_ix)
%         yticklab{f} = num2str(st_tfr.freq(tick_ix(f)),'%.1f');
%     end
    
    % Find max beta points
    beta_pks = zeros([numel(reg_lab) 2]);
    beta_pk_f_ix = zeros([numel(reg_lab) 2]);
    beta_pk_t_ix = zeros([numel(reg_lab) 2]);
    minmax_ix = zeros(size(reg_lab));
    for reg_ix = 1:numel(reg_lab)
        for f = 1:numel(fois)
            if min(beta_mat(reg_ix,f,:)) < beta_pks(reg_ix,1)
                [~, beta_pk_t_ix(reg_ix,1)] = min(beta_mat(reg_ix,f,:));
                beta_pk_f_ix(reg_ix,1) = f;
                beta_pks(reg_ix,1) = beta_mat(reg_ix,beta_pk_f_ix(reg_ix,1),beta_pk_t_ix(reg_ix,1));
            end
            if max(beta_mat(reg_ix,f,:)) > beta_pks(reg_ix,2)
                [~, beta_pk_t_ix(reg_ix,2)] = max(beta_mat(reg_ix,f,:));
                beta_pk_f_ix(reg_ix,2) = f;
                beta_pks(reg_ix,2) = beta_mat(reg_ix,beta_pk_f_ix(reg_ix,2),beta_pk_t_ix(reg_ix,2));
            end
        end
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
    %cfgplt = []; cfgplt.zlim = clim;
    axes = gobjects([numel(reg_lab)+1 1]);
    [num_rc,~] = fn_num_subplots(numel(reg_lab)+1);
    for reg_ix = 1:numel(reg_lab)
        subplot(num_rc(1),num_rc(2),reg_ix);
        axes(reg_ix) = gca; hold on;
        
        % Plot Matrix
        im = imagesc(st_time_vec, fois, squeeze(beta_mat(reg_ix,:,:)),clim);
        im.AlphaData = squeeze(sig_mask(reg_ix,:,:));
        
        % Plot max/min point
        scatter(st_time_vec(beta_pk_t_ix(reg_ix,minmax_ix(reg_ix))),...
            fois(beta_pk_f_ix(reg_ix,minmax_ix(reg_ix))), 100,'r','*');
        
        % Plot Properties
        set(axes(reg_ix),'YDir','normal');
%         set(axes(reg_ix),'YTick',tick_ix);
%         set(axes(reg_ix),'YTickLabels',yticklab);
        set(axes(reg_ix),'YLim',[min(fois) max(fois)]);
        set(axes(reg_ix),'XLim',[min(st_time_vec) max(st_time_vec)]);
        set(axes(reg_ix),'XTick',[min(st_time_vec):0.1:max(st_time_vec)]);
        %ft_singleplotTFR(cfgplt, tfr_avg{cond_ix});
        title([st_tfr.label{ch_ix} ': ' reg_lab{reg_ix} ' Beta']);
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        colorbar('northoutside');
        set(axes(reg_ix),'FontSize',16);
    end
    
    % R2 Plot
    subplot(num_rc(1),num_rc(2),numel(reg_lab)+1);
    axes(numel(reg_lab)+1) = gca; hold on;
    imagesc(st_time_vec, fois, r2_mat);
    set(axes(numel(reg_lab)+1),'YDir','normal');
%     set(axes(numel(reg_lab)+1),'YTick',tick_ix);
%     set(axes(numel(reg_lab)+1),'YTickLabels',yticklab);
    set(axes(numel(reg_lab)+1),'YLim',[min(fois) max(fois)]);
    set(axes(numel(reg_lab)+1),'XLim',[min(st_time_vec) max(st_time_vec)]);
    set(axes(numel(reg_lab)+1),'XTick',[min(st_time_vec):0.1:max(st_time_vec)]);
    %ft_singleplotTFR(cfgplt, tfr_avg{cond_ix});
    title([st_tfr.label{ch_ix} ': R2 (n=' num2str(numel(SBJs)) ')']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    colorbar('northoutside');
    set(gca,'FontSize',16);
    
    % Report min and max beta points
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
