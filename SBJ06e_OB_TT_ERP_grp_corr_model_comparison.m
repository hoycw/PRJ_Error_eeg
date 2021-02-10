function SBJ06e_OB_TT_ERP_grp_corr_model_comparison(SBJ_id,tt_proc_id,ob_proc_id,stat_id,model_id,varargin)
%% Plot OB-TT ERP correlations versus TT RL model regressors per condition
% COMPUTATIONS:
%   Load OB and TT ERP features
%   Load model regressors and average within condition
%   Compute correlations per condition
%   Plot correlations against model regressors
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   tt_proc_id [str] - ID of target time preprocessing pipeline
%   ob_proc_id [str] - ID of oddball preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%       st.model   = feat_id for OB ERP features
%       st.measure = feat_id for TT ERP features
% OUTPUTS:
%   saves scatter plot figure

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
if ~strcmp(st.an_style,'reg'); error('This script is for regression'); end

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
[reg_lab, reg_names, ~, ~]     = fn_regressor_label_styles(model_strs{1});
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(st.stat_cond);

%% Load Oddball ERP Features
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'SBJs');
if numel(SBJs)~=numel(tmp.SBJs) || ~all(strcmp(SBJs,tmp.SBJs)); error('SBJ mismatch!'); end
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'ft_amp','ft_times');

% Z-score feature predictors
ob_amp = ft_amp;

%% Load Target Time ERP Features
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.measure '_' tt_proc_id '.mat'],'SBJs');
if numel(SBJs)~=numel(tmp.SBJs) || ~all(strcmp(SBJs,tmp.SBJs)); error('SBJ mismatch!'); end

% Not doing anything with 'erp_times' right now...
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.measure '_' tt_proc_id '.mat'],...
    'erp_amp','miss_erps');

tt_amp = erp_amp;

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

grp_model = squeeze(mean(model,2));

%% Run Correlation
cond_corr = nan([numel(ft.name) numel(cond_lab)]);
cond_pval = nan([numel(ft.name) numel(cond_lab)]);
for ft_ix = 1:numel(ft.name)
    for cond_ix = 1:numel(cond_lab)
        % Compute correlation
        if any(miss_erps(:,cond_ix))
            fprintf(2,'\tWARNING: %d missing values for %s in %s!\n',...
                sum(miss_erps(:,cond_ix)),st_ft.measure,cond_names{cond_ix});
        end
        [tmp_r,tmp_p] = corrcoef(tt_amp(:,cond_ix),ob_amp(:,ft_ix),'Rows','complete');
        cond_corr(ft_ix,cond_ix) = tmp_r(1,2);
        cond_pval(ft_ix,cond_ix) = tmp_p(1,2);
    end
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(cond_pval,[size(cond_pval,1)*size(cond_pval,2) 1]));
    qvals = reshape(qvals,[size(cond_pval,1) size(cond_pval,2)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

%% Plot Regression results
for ft_ix = 1:numel(ft.name)
    fig_name = [SBJ_id '_' stat_id '_' ft.name{ft_ix} '_' model_id '_amp_compare'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    [n_rowcol,~] = fn_num_subplots(numel(reg_lab));
    for reg_ix = 1:numel(reg_lab)
        subplot(n_rowcol(1),n_rowcol(2),reg_ix); hold on;
                
        % Plot features
        cond_leg  = cell(size(cond_lab));
        cond_scat = gobjects(size(cond_lab));
        for cond_ix = 1:numel(cond_lab)
            cond_leg{cond_ix} = [cond_lab{cond_ix} ' p=' ...
                num2str(cond_pval(ft_ix,cond_ix),'%.3f') ...
                '; q=' num2str(qvals(ft_ix,cond_ix),'%.3f')];
            
            cond_scat(cond_ix) = scatter(grp_model(cond_ix,reg_ix),cond_corr(ft_ix,cond_ix));
%                 mrkr, 'MarkerEdgeColor', cond_colors{cond_ix}, 'MarkerEdgeColor', cond_colors{cond_ix});
            if cond_pval(ft_ix,cond_ix) < 0.01
                cond_scat(cond_ix).Marker = 'p';
                cond_scat(cond_ix).MarkerEdgeColor = cond_colors{cond_ix};
                cond_scat(cond_ix).MarkerFaceColor = cond_colors{cond_ix};
            elseif cond_pval(ft_ix,cond_ix) < 0.05
                cond_scat(cond_ix).Marker = 'd';
                cond_scat(cond_ix).MarkerEdgeColor = cond_colors{cond_ix};
                cond_scat(cond_ix).MarkerFaceColor = cond_colors{cond_ix};
            elseif cond_pval(ft_ix,cond_ix) < 0.1
                cond_scat(cond_ix).Marker = 'o';
                cond_scat(cond_ix).MarkerEdgeColor = cond_colors{cond_ix};
            else
                cond_scat(cond_ix).Marker = '.';
                cond_scat(cond_ix).MarkerEdgeColor = cond_colors{cond_ix};
            end
        end
        
        % Plot parameters
        set(gca,'FontSize',16);
        title([[ft.name{ft_ix} '-' st.measure] ': ' reg_names{reg_ix}]);
        xlabel(reg_lab{reg_ix});
        ylabel('Correlation (r)');
        legend(cond_scat,cond_leg,'Location','best');
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
end

end
