function SBJ06b_OB_ERP_feature_corr(SBJ_id,proc_id,feat_id,varargin)
%% Check correlations between amplitude and latency of condition-averaged ERP peaks in oddball (OB) task
% COMPUTATIONS:
%   Load data and compute correlations across features
%   Plot (scatter) amplitude correlations and latency correlations
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of oddball preprocessing pipeline
%   feat_id [str] - ID of the feature extraction parameters to use
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

%% Load Data 
if ~contains(proc_id,'odd'); error('proc_id must be for oddball task!'); end
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' feat_id '_vars.m'];
eval(feat_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(ft.grp_id);

%% Load Features
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' feat_id '_' proc_id '.mat']);

%% Compare Amplitude Features
fig_name = ['ERPs_' feat_id '_' SBJ_id '_amp_corr'];
figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.7],'Visible',fig_vis);
[n_rowcol,~] = fn_num_subplots(nchoosek(numel(ft.name),2));

pair_ix = 0;
for ft_ix1 = 1:numel(ft.name)
    for ft_ix2 = ft_ix1+1:numel(ft.name)
        pair_ix = pair_ix + 1;
        subplot(n_rowcol(1),n_rowcol(2),pair_ix); hold on;
        
        % Compute correlation
        [r,p] = corrcoef(ft_amp(ft_ix1,:),ft_amp(ft_ix2,:));
        r = r(1,2); p = p(1,2);
        
        % Plot features
        scatter(ft_amp(ft_ix1,:),ft_amp(ft_ix2,:), 'o', 'k');
        
        % Plot linear fit
        coeff = polyfit(ft_amp(ft_ix1,:),ft_amp(ft_ix2,:),1);
        xbounds = get(gca,'XLim');
        xdat = [xbounds(1)+1 xbounds(2)-1];
        ydat = coeff(1)*xdat + coeff(2);
        line(xdat,ydat);
        
        
        % Plot parameters
        title(['r=' num2str(r,'%.3f') '; p=' num2str(p,'%.3f')]);
        xlabel([ft.name{ft_ix1} '(' ft.chan{ft_ix1} ', ' ...
            ft.cond{ft_ix1} ') Amp (uV)']);
        ylabel([ft.name{ft_ix2} '(' ft.chan{ft_ix2} ', ' ...
            ft.cond{ft_ix2} ') Amp (uV)']);
        set(gca,'FontSize',16);
    end
end

% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' ft.an_id '/' feat_id '/'];
    if ~exist(fig_dir,'dir') && save_fig
        mkdir(fig_dir);
    end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Compare Latency Features
if ~strcmp(ft.measure,'grpMW')
    fig_name = ['ERPs_' feat_id '_' SBJ_id '_lat_corr'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.7],'Visible',fig_vis);
    [n_rowcol,~] = fn_num_subplots(nchoosek(numel(ft.name),2));
    
    pair_ix = 0;
    for ft_ix1 = 1:numel(ft.name)
        for ft_ix2 = ft_ix1+1:numel(ft.name)
            pair_ix = pair_ix + 1;
            subplot(n_rowcol(1),n_rowcol(2),pair_ix); hold on;
            
            % Compute correlation
            [r,p] = corrcoef(ft_times(ft_ix1,:),ft_times(ft_ix2,:));
            r = r(1,2); p = p(1,2);
            
            % Plot features
            scatter(ft_times(ft_ix1,:),ft_times(ft_ix2,:), 'o', 'k');
            
            % Plot linear fit
            coeff = polyfit(ft_times(ft_ix1,:),ft_times(ft_ix2,:),1);
            xbounds = get(gca,'XLim');
            xdat = [xbounds(1)+0.01 xbounds(2)-0.01];
            ydat = coeff(1)*xdat + coeff(2);
            line(xdat,ydat);
            
            
            % Plot parameters
            title(['r=' num2str(r,'%.3f') '; p=' num2str(p,'%.3f')]);
            xlabel([ft.name{ft_ix1} '(' ft.chan{ft_ix1} ', ' ...
                ft.cond{ft_ix1} ') Lat (s)']);
            ylabel([ft.name{ft_ix2} '(' ft.chan{ft_ix2} ', ' ...
                ft.cond{ft_ix2} ') Lat (s)']);
            set(gca,'FontSize',16);
        end
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Compute variance inflation factors
amp_vifs = fn_variance_inflation_factor(zscore(ft_amp'));
lat_vifs = fn_variance_inflation_factor(zscore(ft_times'));

for ft_ix = 1:numel(ft.name)
    if amp_vifs(ft_ix)>10
        fprintf(2,'%s: %s (%s, %s) amp;lat VIF = %.2f; %.2f\n',SBJ_id,ft.measure,...
            ft.name{ft_ix},ft.cond{ft_ix},amp_vifs(ft_ix),lat_vifs(ft_ix));
    else
        fprintf('%s: %s (%s, %s) amp;lat VIF = %.2f; %.2f\n',SBJ_id,ft.measure,...
            ft.name{ft_ix},ft.cond{ft_ix},amp_vifs(ft_ix),lat_vifs(ft_ix));
    end
end

end
