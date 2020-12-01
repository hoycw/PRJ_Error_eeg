function SBJ06c_OB_ERP_TT_feature_corr(SBJ_id,ob_proc_id,tt_proc_id,an_id,ob_feat_id,stat_id,varargin)
%% Compute correlations between oddball and target time ERP features
% COMPUTATIONS:
%   Load data and compute correlations across features
%   Plot (scatter) amplitude correlations and latency correlations
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   ob_proc_id [str] - ID of oddball preprocessing pipeline
%   tt_proc_id [str] - ID of target time preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   ob_feat_id [str] - ID of the OB feature extraction parameters to use
%   tt_stat_id [str] - ID of the TT stat parameters to use
%   varargin:
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   save figures of correlations between oddball ERP features

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
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' ob_feat_id '/'];
if ~exist(fig_dir,'dir') && save_fig
    mkdir(fig_dir);
end

%% Load Data 
if ~strcmp(ob_proc_id,'odd_full_ft'); error('why not odd_full_ft for oddball task?'); end
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' ob_proc_id '_vars.m'];
eval(proc_vars_cmd);
if ~strcmp(ob_proc_id,'eeg_full_ft'); error('why not eeg_full_ft for target time task?'); end
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' tt_proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' ob_feat_id '_vars.m'];
eval(feat_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(ft.conditions);

%% Load Features
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' ob_feat_id '_' an_id '.mat']);

%% Compare Amplitude Features
fig_name = ['ERPs_' ob_feat_id '_' SBJ_id '_amp_corr'];
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
fig_name = ['ERPs_' ob_feat_id '_' SBJ_id '_lat_corr'];
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
