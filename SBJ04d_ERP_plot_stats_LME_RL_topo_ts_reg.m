function SBJ04d_ERP_plot_stats_LME_RL_topo_ts_reg(SBJ_id,an_id,stat_ids,plt_id,save_fig,varargin)
%% Plots group RL beta topographies per regressor across multiple stats run
% at different time points
%   IMPORTANT: Recomputes FDR correction over all analyses (times), regressors, and electrodes
%   Rows are regressors, columns are time points
%   Prints statistics: max beta and q value per regressor across time points
%       and electrodes, along with peak electrode and time
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   an_id [str] - ID of the analysis parameters to use
%   stat_ids [cell array] - string IDs of the stats parameters to compare
%   plt_id [str] - ID of the plotting parameters to use
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

%% Analysis and Plotting Parameters
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
if ~strcmp(an.ROI{1},'all'); error('run this only for topos with all channels!'); end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Load stat parameters and check compatibility
sts        = cell(size(stat_ids));
model_labs = cell(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    % Load stat_id
    stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_ids{st_ix} '_vars.m'];
    eval(stat_vars_cmd);
    sts{st_ix} = st;
    model_labs{st_ix} = sts{st_ix}.model_lab;
    if ~strcmp(st.measure,'mean'); error('this script is for mean window only!');end
    
    % Ensure all stat params are equal except peak center
    if st_ix>1
        fnames = fieldnames(sts{st_ix});
        for f_ix = 1:numel(fnames)
            if ~any(strcmp(fnames{f_ix},{'pk_center','pk_reg_id','pk_an_id'}))%,'stat_cond'
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

% Select Conditions of Interest
[reg_lab, reg_names, ~, ~]  = fn_regressor_label_styles(sts{1}.model_lab);

%% Load Stats
clim     = zeros([numel(reg_lab) 2]);
betas    = nan([numel(stat_ids) numel(reg_lab) 64]);
pvals    = nan([numel(stat_ids) numel(reg_lab) 64]);
r2       = zeros([numel(stat_ids) 64]);
pk_times = zeros(size(stat_ids));
for st_ix = 1:numel(stat_ids)
    load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_ids{st_ix} '_' an_id '.mat'],'lme','qvals','ch_list','reg_pk_time');
    if numel(ch_list)<64; error('Cannot plot topo without full cap!'); end
    
    % Get beta values, p values, color limits, and peak times
    pk_times(st_ix) = reg_pk_time;
    for ch_ix = 1:numel(ch_list)
        for reg_ix = 1:numel(reg_lab)
            betas(st_ix,reg_ix,ch_ix) = lme{ch_ix}.Coefficients.Estimate(reg_ix+1);
            clim(reg_ix,:) = [min([clim(reg_ix,1) betas(st_ix,reg_ix,ch_ix)]) max([clim(reg_ix,2) betas(st_ix,reg_ix,ch_ix)])];
        end
        pvals(st_ix,:,ch_ix) = lme{ch_ix}.Coefficients.pValue(2:end);
        r2(st_ix,ch_ix) = lme{ch_ix}.Rsquared.Adjusted;
    end
    
    clear lme qvals reg_pk_time
end

% Recompute significance with multiple comparisons over all regressors and windows
if strcmp(sts{1}.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2)*size(pvals,3) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2) size(pvals,3)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end
sig_ch = qvals<=sts{1}.alpha;

%% Create dummy dataset for plotting
% This will initialize electrode labels, neighbor structure, etc.
load([root_dir 'PRJ_Error_eeg/data/',SBJs{1},'/04_proc/',SBJs{1},'_',an_id,'.mat'],'roi');
topo = {};
topo.label  = roi.label;
topo.dimord = 'chan_time';

%% Plot Results
% Create plot
fig_name = [SBJ_id '_topo_ts_' an_id];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);

% Set up defatul plotting params
axes = gobjects([numel(stat_ids) numel(reg_lab)+1]);
cfgp = [];
cfgp.layout          = 'biosemi64.lay';
cfgp.comment         = 'no';
cfgp.highlight       = 'on';
cfgp.highlightsymbol = '*';
cfgp.parameter       = 'avg';
% cfgp.maskparameter = 'mask';
cfgp.zlim            = [-max(abs(clim(:))) max(abs(clim(:)))];
cfgp.colormap        = 'jet';

% Plot topo for each regressor and time point
plot_ix = 0;
for reg_ix = 1:numel(reg_lab)
    for st_ix = 1:numel(stat_ids)
        topo.time   = pk_times(st_ix);
        cfgp.xlim   = [pk_times(st_ix) pk_times(st_ix)];
        plot_ix = plot_ix + 1;
        subplot(numel(reg_lab)+1,numel(stat_ids),plot_ix);
        axes(st_ix,reg_ix) = gca; hold on;
        
        % Only plot color bar on last column
        if st_ix==numel(stat_ids)
            cfgp.colorbar = 'yes';
        else
            cfgp.colorbar = 'no';
        end
        
        % Plot Beta Topos
        topo.avg  = squeeze(betas(st_ix,reg_ix,:));
        cfgp.highlightchannel = find(sig_ch(st_ix,reg_ix,:));
        ft_topoplotER(cfgp, topo);
        title([reg_names{reg_ix} ': ' num2str(pk_times(st_ix),'%.3f') ' s']);
        axis tight
    end
end

% Plot R2 as bottom row
cfgp.zlim = [0 max(r2(:))];
cfgp.highlight = 'no';
cfgp.highlightchannel = [];
for st_ix = 1:numel(stat_ids)
    plot_ix = plot_ix + 1;
    subplot(numel(reg_lab)+1,numel(stat_ids),plot_ix);
    axes(st_ix,numel(reg_ix)+1) = gca; hold on;
    
    % Only plot color bar on last column
    if st_ix==numel(stat_ids)
        cfgp.colorbar = 'yes';
    else
        cfgp.colorbar = 'no';
    end
    
    % Plot R2 Topos
    topo.avg  = r2(st_ix,:)';
    ft_topoplotER(cfgp, topo);
    title(['Adjusted R2: '  num2str(pk_times(st_ix),'%.3f') ' s']);
    axis tight    
end

%% Report peak beta, q value, window, and elec per regressor
for reg_ix = 1:numel(reg_lab)
    max_beta = 0; max_ch_ix = 0; max_t_ix = 0;
    for st_ix = 1:numel(stat_ids)
        if max(abs(betas(st_ix,reg_ix,:))) > abs(max_beta)
            max_tmp = max(betas(st_ix,reg_ix,:));
            min_tmp = min(betas(st_ix,reg_ix,:));
            if abs(max_tmp) > abs(min_tmp)
                [max_beta, max_ch_ix] = max(betas(st_ix,reg_ix,:));
            else
                [max_beta, max_ch_ix] = min(betas(st_ix,reg_ix,:));
            end
            max_t_ix = st_ix;
        end
    end
    fprintf('%s max beta = %.03f at %.03f, %s; p = %.20f\n',reg_lab{reg_ix},max_beta,...
        pk_times(max_t_ix),roi.label{max_ch_ix},qvals(max_t_ix,reg_ix,max_ch_ix));
end

%% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' strjoin(stat_ids,'-') '/' plt_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    % Ensure vector graphics if saving
    if any(strcmp(fig_ftype,{'svg','eps'}))
        set(gcf, 'Renderer', 'painters');
    end
    saveas(gcf,fig_fname);
end

end
