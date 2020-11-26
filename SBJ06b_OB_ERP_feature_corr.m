function SBJ06b_OB_ERP_feature_corr(SBJ_id,proc_id,an_id,feat_id,varargin)
%% Check correlations between amplitude and latency of condition-averaged ERP peaks in oddball (OB) task
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load and average single-trial ERP and design matrix (model regressors, SBJ factor) within condition
%       Optional: z-score model regressors within SBJ
%   Compute and plot group concatenated model (design matrix) and correlations
%   Find positive and negative peaks:
%       Must find positive and negative peaks in their respective windows
%       Positive peak must precede negative peak
%   Compute amplitude differences:
%       For multiple peak pairs, selects largest amplitude difference
%       Rejects pairs with amplitude difference in wrong direction (negative > positive)
%   Plot peak detection errors
%       Optional: Plot ERPs and detected peaks
%   Run linear mixed effects model per time point or per electrode
%   Correct for multiple comparisons (FDR for regressors)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   feat_id [str] - ID of the stats parameters to use
%   varargin:
%       plot_feat [0/1] - binary flag for plotting each ERP with peaks overlay
%           default: 1
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)
%   pk_amp [float] - amplitudes of [positive, negative] peaks
%   pk_times [float] - times (in sec) of [positive, negative] peaks

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

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
if ~exist('fig_vis','var');    fig_vis = 'on'; end
if ~exist('fig_ftype','var');  fig_ftype = 'png'; end
if ~exist('save_fig','var');   save_fig = 1; end
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' feat_id '/'];
if ~exist(fig_dir,'dir') && save_fig
    mkdir(fig_dir);
end

%% Load Data 
if ~contains(proc_id,'odd'); error('proc_id must be for oddball task!'); end
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' feat_id '_vars.m'];
eval(feat_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(ft.conditions);

%% Load Features
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' feat_id '_' an_id '.mat']);

%% Compare Amplitude Features
fig_name = ['ERPs_' feat_id '_' SBJ_id '_amp_corr'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.7 1],'Visible',fig_vis);
[n_rowcol,~] = fn_num_subplots(nchoosek(numel(ft.feat_name),2));

pair_ix = 0;
for feat_ix1 = 1:numel(ft.feat_name)
    for feat_ix2 = feat_ix1+1:numel(ft.feat_name)
        pair_ix = pair_ix + 1;
        subplot(n_rowcol(1),n_rowcol(2),pair_ix); hold on;
        
        % Compute correlation
        [r,p] = corrcoef(pk_amp(feat_ix1,:),pk_amp(feat_ix2,:));
        r = r(1,2); p = p(1,2);
        
        % Plot features
        scatter(pk_amp(feat_ix1,:),pk_amp(feat_ix2,:), 'o', 'k');
        
        % Plot linear fit
        coeff = polyfit(pk_amp(feat_ix1,:),pk_amp(feat_ix2,:),1);
        xbounds = get(gca,'XLim');
        xdat = [xbounds(1)+1 xbounds(2)-1];
        ydat = coeff(1)*xdat + coeff(2);
        line(xdat,ydat);
        
        
        % Plot parameters
        title(['r=' num2str(r,'%.3f') '; p=' num2str(p,'%.3f')]);
        xlabel([ft.feat_name{feat_ix1} '(' ft.pk_chan{feat_ix1} ', ' ...
            ft.feat_cond{feat_ix1} ') Amp (uV)']);
        ylabel([ft.feat_name{feat_ix2} '(' ft.pk_chan{feat_ix2} ', ' ...
            ft.feat_cond{feat_ix2} ') Amp (uV)']);
        set(gca,'FontSize',16);
    end
end

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compare Latency Features
fig_name = ['ERPs_' feat_id '_' SBJ_id '_lat_corr'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.7 1],'Visible',fig_vis);
[n_rowcol,~] = fn_num_subplots(nchoosek(numel(ft.feat_name),2));

pair_ix = 0;
for feat_ix1 = 1:numel(ft.feat_name)
    for feat_ix2 = feat_ix1+1:numel(ft.feat_name)
        pair_ix = pair_ix + 1;
        subplot(n_rowcol(1),n_rowcol(2),pair_ix); hold on;
        
        % Compute correlation
        [r,p] = corrcoef(pk_times(feat_ix1,:),pk_times(feat_ix2,:));
        r = r(1,2); p = p(1,2);
        
        % Plot features
        scatter(pk_times(feat_ix1,:),pk_times(feat_ix2,:), 'o', 'k');
        
        % Plot linear fit
        coeff = polyfit(pk_times(feat_ix1,:),pk_times(feat_ix2,:),1);
        xbounds = get(gca,'XLim');
        xdat = [xbounds(1)+0.01 xbounds(2)-0.01];
        ydat = coeff(1)*xdat + coeff(2);
        line(xdat,ydat);
        
        
        % Plot parameters
        title(['r=' num2str(r,'%.3f') '; p=' num2str(p,'%.3f')]);
        xlabel([ft.feat_name{feat_ix1} '(' ft.pk_chan{feat_ix1} ', ' ...
            ft.feat_cond{feat_ix1} ') Lat (s)']);
        ylabel([ft.feat_name{feat_ix2} '(' ft.pk_chan{feat_ix2} ', ' ...
            ft.feat_cond{feat_ix2} ') Lat (s)']);
        set(gca,'FontSize',16);
    end
end

% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compute variance inflation factors
%   This is useless: "Warning: Matrix is close to singular or badly scaled.
%   Results may be inaccurate."
% amp_vifs = fn_variance_inflation_factor(zscore(pk_amp));
% lat_vifs = fn_variance_inflation_factor(zscore(pk_times));

end
