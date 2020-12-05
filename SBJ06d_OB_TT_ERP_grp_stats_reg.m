function SBJ06d_OB_TT_ERP_grp_stats_reg(SBJ_id,tt_proc_id,ob_proc_id,stat_id,varargin)
%% Run multiple regression using OB ERP features to predict TT ERP features
%   Oddball ERPs only yield one feature per SBJ, so can only predict single TT condition/average
%   If multiple TT conditions, then predicting each one separately
% COMPUTATIONS:
%   Load OB and TT ERP features
%       Optional: z-score regressors across group
%   Run general linear model with OB features predicting each TT condition
%   Correct for multiple comparisons (FDR for regressors)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   tt_proc_id [str] - ID of target time preprocessing pipeline
%   ob_proc_id [str] - ID of oddball preprocessing pipeline
%   stat_id [str] - ID of the stats parameters to use
%       st.model   = feat_id for OB ERP features
%       st.measure = feat_id for TT ERP features
% OUTPUTS:
%   glm [cell array] - GeneralLinearModel output class, one cell per channel
%   qvals [float array] - [n_regressors, n_chan/n_time] p values adjusted for multiple comparisons 
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)

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
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.an_style,'reg'); error('This script is for regression'); end

% TT Feature Parameters
stat_feat_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' st.measure '_vars.m'];
eval(stat_feat_cmd);
st_ft = ft;
if ~any(strcmp(st_ft.feat_lab,{'FRN','P3'})); error('This script is only ready for FRN and P3!'); end
if ~any(strcmp(st_ft.grp_id,{'All','DifFB','Pos','Neg'})); error('stat feat should be TT conditions!'); end

% OB Feature Parameters
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' st.model_lab '_vars.m'];
eval(feat_vars_cmd);
if ~any(strcmp(ft.grp_id,{'rare','Odd','Tar'})); error('Features should be oddball conditions!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);

%% Load Oddball ERP Features
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'SBJs');
if numel(SBJs)~=numel(tmp.SBJs) || ~all(strcmp(SBJs,tmp.SBJs)); error('SBJ mismatch!'); end
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.model_lab '_' ob_proc_id '.mat'],'ft_amp','ft_times');

% Z-score feature predictors
model = ft_amp;
if st.z_reg
    model = zscore(model);
end

%% Load Target Time ERP Features
tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.measure '_' tt_proc_id '.mat'],'SBJs');
if numel(SBJs)~=numel(tmp.SBJs) || ~all(strcmp(SBJs,tmp.SBJs)); error('SBJ mismatch!'); end
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' st.measure '_' tt_proc_id '.mat'],...
    'erp_amp','erp_times','miss_erps');
data = erp_amp;

%% Run Linear Multiple Regression
fprintf('========================== Running Stats ==========================\n');
% Build Model Table
tbl = table;
for ft_ix = 1:numel(ft.name)
    tbl.(ft.name{ft_ix}) = model(:,ft_ix);
end

% Create Model Formula
reg_formula = strjoin(ft.name,' + ');
formula = ['ERP ~ ' reg_formula];

% Run Model: Single time point (mean in window) over channels
glm = cell(size(cond_lab));
pvals = nan([numel(ft.name) numel(cond_lab)]);
for cond_ix = 1:numel(cond_lab)
    tbl.ERP = data(:,cond_ix);
    glm{cond_ix} = fitlme(tbl,formula);
    pvals(:,cond_ix) = glm{cond_ix}.Coefficients.pValue(2:end); % Skip intercept
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

%% Plot Regression results
for ft_ix = 1:numel(ft.name)
    fig_name = [SBJ_id '_' stat_id '_' ft.name{ft_ix} '_amp_fits'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    [n_rowcol,~] = fn_num_subplots(numel(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        subplot(n_rowcol(1),n_rowcol(2),cond_ix); hold on;
        
        % Compute correlation
        if any(miss_erps(:,cond_ix))
            fprintf(2,'\tWARNING: %d missing values for %s in %s!\n',...
                sum(miss_erps(:,cond_ix)),st_ft.measure,cond_lab{cond_ix});
        end
        [r,p] = corrcoef(model(:,ft_ix),data(:,cond_ix),'Rows','complete');
        r = r(1,2); p = p(1,2);
        
        % Plot features
        scatter(model(:,ft_ix),data(:,cond_ix), 'o', 'k');
        
        % Plot linear fit
        coeff = polyfit(model(:,ft_ix),data(:,cond_ix),1);
        xbounds = get(gca,'XLim');
        xfudge = (xbounds(2)-xbounds(1))*0.1;
        xdat = [xbounds(1)+xfudge xbounds(2)-xfudge];
        ydat = coeff(1)*xdat + coeff(2);
        simple_fit = line(xdat,ydat);
        
        % Plot parameters
        set(gca,'FontSize',16);
        legend(simple_fit,['r=' num2str(r,'%.2f') '; p=' num2str(p,'%.3f')],...
            'Location','best');
        title([cond_lab{cond_ix} ': beta=' num2str(glm{cond_ix}.Coefficients.Estimate(ft_ix+1)) ...
            '; q=' num2str(qvals(ft_ix,cond_ix),'%.3f')]);
        xlabel([ft.name{ft_ix} '(' ft.chan{ft_ix} ', ' ft.cond{ft_ix} ') Amp (uV)']);
        ylabel([st_ft.name{1} ' ' st_ft.measure ' Amp (uv)']);
    end
    
    % Save figure
    if save_fig
        fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' st_ft.an_id '/' stat_id '/'];
        if ~exist(fig_dir,'dir') && save_fig
            mkdir(fig_dir);
        end
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Save Results
stat_out_fname = [root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','glm','qvals','SBJs');

end
