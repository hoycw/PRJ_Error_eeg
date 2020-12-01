%% Extract and compare features of Oddball and Target Time ERPs for Sequential PE Revision
% Started 11/25/20 by Colin W Hoy
% Final model_id = 'ERPEsL_all' (all = 'DifFB', which includes all conditions)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%%
addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% General parameters
SBJ_id = 'goodEEG1';%'goodEEG';%'goodEEG2';%
SBJs = fn_load_SBJ_list(SBJ_id);

%% Oddball ERP feature extraction
proc_id  = 'odd_full_ft';
feat_ids = {'N2bP3a_sbjPk','N2bP3a_sbjMW1','N2bP3a_grpMW1'};
%   an_id is specified in the feat struct (always 'ERP_all_S2t1_dm2t0_fl05t20')

for ft_ix = 1:numel(feat_ids)
    % Extract Oddball ERP features
    SBJ06a_OB_ERP_save_features(SBJ_id,proc_id,feat_ids{ft_ix});
    
    % Plot correlations between ERP features before using as model predictors
    SBJ06b_OB_ERP_feature_corr(SBJ_id,proc_id,feat_ids{ft_ix});
end

%% Oddball vs. Target Time ERP Comparison
proc_id  = 'eeg_full_ft';
an_id    = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_ids = {'OB_N2P3_DifFB_reg_erpmn1FRN'};

SBJ06c_OB_ERP_TT_grp_stats_LME_mean_window(SBJ_id,proc_id,an_id,stat_id);


plt_id    = 'bar_sigStar';
null_id   = 'SBJonly_all_lme_mn1FRN';

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

for st_ix = 1:numel(stat_ids)
    if ~isempty(strfind(stat_ids{st_ix},'erpmn'))
        % Average across ERPs (e.g., stat_id = 'ERPEsL_all_lme_erpmn1FRN')
        SBJ04c_ERP_grp_stats_LME_mean_window(SBJ_id,proc_id,an_id,stat_ids{st_ix});
    else
        % Average across single trials (e.g., stat_id = 'ERPEsL_all_lme_mn1FRN')
        %   Not used because literature typically averages over ERPs, not
        %   single trials, hence ERP mean
%         SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_id,stat_ids{st_ix});
    end

    % Sup. Fig. 1C: Plot mean window betas
    SBJ04d_ERP_plot_stats_LME_mean_betas(SBJ_id,proc_id,an_id,stat_ids{st_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

    % Plot mean window betas with data as bar plot
%     SBJ04d_ERP_plot_stats_LME_mean_betas(SBJ_id,proc_id,an_id,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_data',1);
end

