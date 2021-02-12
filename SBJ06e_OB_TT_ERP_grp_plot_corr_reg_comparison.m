function SBJ06e_OB_TT_ERP_grp_plot_corr_reg_comparison(SBJ_id,tt_proc_id,stat_id,model_id,varargin)
%% Plot OB-TT ERP correlations versus TT RL model regressors per condition
% COMPUTATIONS:
%   Load OB-TT ERP correlation results from SBJ06d
%   Load model regressors and average within condition
%   Plot OB-TT correlations against model regressors (condition subplot)
%   Plot OB-TT correlation matrix (TT conditions, OB features)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   tt_proc_id [str] - ID of target time preprocessing pipeline
%   ob_proc_id [str] - ID of oddball preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%       st.model   = feat_id for OB ERP features
%       st.measure = feat_id for TT ERP features
% OUTPUTS:
%   saves OB-TT corr vs. RL model reg scatter plot figure

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
        elseif strcmp(varargin{v},'save_fig')
            save_fig = varargin{v+1};
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
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/OB_TT_feat/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.an_style,'corr'); error('This script is for OB-TT correlation'); end

% TT Feature Parameters
stat_feat_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' st.measure '_vars.m'];
eval(stat_feat_cmd);
st_ft = ft;
if ~any(strcmp(st_ft.name,{'FRN','P3','sRPE','uRPE','Lik'}))
    error('This script is only ready for FRN, P3, and reg peaks!');
end
if ~any(strcmp(st_ft.grp_id,{'All','DifFB','Pos','Neg'})); error('stat feat should be TT conditions!'); end

% OB Feature Parameters
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' st.model_lab '_vars.m'];
eval(feat_vars_cmd);
if ~any(strcmp(ft.grp_id,{'rare','Odd','Tar'})); error('Features should be oddball conditions!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
model_strs = strsplit(model_id,'_');
[reg_lab, reg_names, reg_colors, ~]     = fn_regressor_label_styles(model_strs{1});
[cond_lab, cond_names, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);
[ob_ft_lab, ob_ft_names, ob_ft_colors, ob_ft_styles, ob_ft_mrkrs] = fn_OB_ERP_label_styles(ft.feat_lab);
if numel(ob_ft_lab)~=numel(ft.name) || ~all(strcmp(ob_ft_lab,ft.name)); error('OB feature mismatch!'); end

%% Load Correlation Stat Output
stat_fname = [root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '.mat'];
tmp = load(stat_fname,'SBJs');
if numel(SBJs)~=numel(tmp.SBJs) || ~all(strcmp(SBJs,tmp.SBJs)); error('SBJ mismatch!'); end

load(stat_fname,'cond_corr','cond_pval','cond_qval');

%% Load Data and Build Model
model = nan([numel(cond_lab) numel(SBJs) numel(reg_lab)]);
for s = 1:numel(SBJs)
    % Load behavior
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' tt_proc_id '_final.mat'],'bhv');
    cond_idx = fn_condition_index(cond_lab, bhv);
    
    % Load model (no z-score)
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
    
    % Average model within condition
    for cond_ix = 1:numel(cond_lab)
        model(cond_ix,s, :) = mean(tmp.model(cond_idx==cond_ix,:),1);
    end
    
    clear bhv tmp cond_idx
end

grp_model_avg = squeeze(nanmean(model,2));
grp_model_std = squeeze(nanstd(model,[],2));
grp_model_sem = grp_model_std./sqrt(size(model,2))';

%% Plot Correlation results against Regressor Values
fig_name = [SBJ_id '_' stat_id '_' model_id '_corr_reg_compare'];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
[n_rowcol,~] = fn_num_subplots(numel(reg_lab));

% Create significance scatter size mapping
sig_sz_min = 25;
sig_sz_max = 100;
sig_sz_step = 100;
sig_szs    = linspace(sig_sz_min,sig_sz_max,sig_sz_step);
sig_sz_map = linspace(0.05,min(cond_qval(:)),sig_sz_step);

for reg_ix = 1:numel(reg_lab)
    subplot(n_rowcol(1),n_rowcol(2),reg_ix);
    ax1 = gca; hold on;
    
    ft_lines = gobjects(size(ft.name));
    ft_leg   = cell(size(ft.name));
    for ft_ix = 1:numel(ft.name)
        % Plot Condition Correlations
        for cond_ix = 1:numel(cond_lab)
            [~,sig_sz_ix] = min(abs(cond_qval(ft_ix,cond_ix)-sig_sz_map));
            if cond_qval(ft_ix,cond_ix)<=0.05; mrkr = ob_ft_mrkrs{ft_ix}; else; mrkr = '.'; end
            scatter(cond_ix,cond_corr(ft_ix,cond_ix),...
                sig_szs(sig_sz_ix), mrkr, 'MarkerEdgeColor', ob_ft_colors{ft_ix});
        end
        ft_lines(ft_ix) = line(1:numel(cond_lab),cond_corr(ft_ix,:),...
            'Color',ob_ft_colors{ft_ix},'LineStyle',ob_ft_styles{ft_ix},'LineWidth',2);
        
%         % Plot linear fit
%         coeff = polyfit(grp_model_avg(:,reg_ix),cond_corr(ft_ix,:)',1);
%         xbounds = get(gca,'XLim');
%         xfudge = (xbounds(2)-xbounds(1))*0.1;
%         xdat = [xbounds(1)+xfudge xbounds(2)-xfudge];
%         ydat = coeff(1)*xdat + coeff(2);
%         ft_lines(ft_ix) = line(xdat,ydat,'Color',ob_ft_colors{ft_ix},'LineWidth',2);
        
        [r,p] = corrcoef(grp_model_avg(:,reg_ix),cond_corr(ft_ix,:)');
        if p(1,2)<0.05; ft_sig_str = '*'; else; ft_sig_str = ''; end
        ft_leg{ft_ix} = [ft_sig_str ob_ft_names{ft_ix} '(p=' num2str(p(1,2),'%.3f') ')'];
    end
    ax1.XLabel.String = 'Conditions';
    ax1.XLim          = [0.5 numel(cond_lab)+0.5];
    ax1.XTick         = 1:numel(cond_lab);
    ax1.XTickLabel    = cond_names;
    ax1.YLim          = [0 1];
    ax1.YLabel.String = 'Correlation (r)';
    ax1.FontSize      = 16;
%     ax1.Color         = [0 0 0];
    
    % Plot Mean Regressor with SEM
    yyaxis right
    ax2 = gca;
    reg_line = errorbar(grp_model_avg(:,reg_ix),grp_model_sem(:,reg_ix),...
                        'Color',reg_colors{reg_ix}, 'LineWidth',3);
    ax2.YLabel.String = reg_names{reg_ix};
    ax2.YColor        = reg_colors{reg_ix};
    ax2.FontSize      = 16;
    
    % Plot parameters
    title([st_ft.chan{1} ' ' st_ft.pk_model_lab ' ' st_ft.pk_reg_id ' Win: ' reg_names{reg_ix}]);
    legend([ft_lines reg_line],[ft_leg reg_lab(reg_ix)],'Location','best');
end

% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' st_ft.an_id '/' stat_id '/' model_id '/'];
    if ~exist(fig_dir,'dir') && save_fig
        mkdir(fig_dir);
    end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

%% Plot Correlation Matrix Summary
fig_name = [SBJ_id '_' stat_id '_corr_mat'];
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
hold on;

% Plot Correlation Matrix
imagesc(cond_corr);
colorbar;
caxis([0 1]);
% colormap(redblue); caxis([-1 1]);

% Plot significance
sig_xfudge = 0.1;
sig_cuts   = [0.05 0.01 0.001];
for ft_ix = 1:numel(ft.name)
    for cond_ix = 1:numel(cond_lab)
        if cond_qval(ft_ix,cond_ix)<=sig_cuts(1)
            if cond_qval(ft_ix,cond_ix)<=sig_cuts(3)
                x_pos = [-sig_xfudge 0 sig_xfudge];
            elseif cond_qval(ft_ix,cond_ix)<=sig_cuts(2)
                x_pos = [-sig_xfudge sig_xfudge]/2;
            else
                x_pos = 0;
            end
            scatter(x_pos+cond_ix,repmat(ft_ix,size(x_pos)),100,'k','*');
        end
    end
end

% Plot Properties
set(gca,'XLim',[0.5 numel(cond_lab)+0.5]);
set(gca,'XTickLabel',cond_names);
set(gca,'YLim',[0.5 numel(ft.name)+0.5]);
set(gca,'YTick',1:numel(ft.name));
set(gca,'YTickLabel',ob_ft_names);
title([st_ft.chan{1} ': ' st_ft.pk_model_lab ' ' st_ft.feat_lab ' (pFDR = 3*<' num2str(sig_cuts(3),'%.03f') '; 2*<' ...
    num2str(sig_cuts(2),'%.02f') '; *<' num2str(sig_cuts(1),'%.02f') ')'], 'interpreter', 'none');
set(gca,'FontSize',16);

% Save figure
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' st_ft.an_id '/' stat_id '/' model_id '/'];
    if ~exist(fig_dir,'dir') && save_fig
        mkdir(fig_dir);
    end
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    saveas(gcf,fig_fname);
end

end
