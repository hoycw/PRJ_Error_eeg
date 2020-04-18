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
SBJ_id = 'good1';
SBJs = load_SBJ_file(SBJ_id);

%% Save Group ERP Peak Time
conditions = 'All';
an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
proc_id    = 'eeg_full_ft';

SBJ03c_ERP_save_grp_ERP_cond(SBJ_id,conditions,proc_id,an_id);

%% Single SBJ RL Model
proc_id  = 'eeg_full_ft';
stat_ids = {'RV_all_lme_mn1FRN','RVL_all_lme_mn1FRN','RVM_all_lme_mn1FRN','RVLM_all_lme_mn1FRN',...
    'sRPE_all_lme_mn1FRN','uRPE_all_lme_mn1FRN','RPE_all_lme_mn1FRN','RPEOL_all_lme_mn1FRN'};

for s = 1:numel(SBJs)
    for st_ix = 1:numel(stat_ids)
        SBJ04a_RL_model(SBJs{s},proc_id,stat_ids{st_ix});
    end
    close all;
end

%% ERP: Linear Mixed Effects Model (Over Time)
% Main RL Model
an_id     = 'ERP_Fz_F2t1_dm2t0_fl05t20';
% stat_ids  = {'RV_all_lme_mn1FRN','RVL_all_lme_mn1FRN','RVM_all_lme_mn1FRN','RVLM_all_lme_mn1FRN'};
% stat_ids  = {'sRPE_all_lme_mn1FRN','uRPE_all_lme_mn1FRN','RPE_all_lme_mn1FRN','RPEOL_all_lme_mn1FRN'};
stat_ids  = {'RV_all_lme_mn1FRN','RVL_all_lme_mn1FRN','RVM_all_lme_mn1FRN','RVLM_all_lme_mn1FRN',...
    'sRPE_all_lme_mn1FRN','uRPE_all_lme_mn1FRN','RPE_all_lme_mn1FRN','RPEOL_all_lme_mn1FRN'};
plt_id    = 'bar_sigStar';
null_id   = 'SBJonly_all_lme_mn1FRN';

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for st_ix = 1:numel(stat_ids)
%     SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_id,stat_ids{st_ix});
%     SBJ04d_ERP_plot_stats_LME_betas(SBJ_id,an_id,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end
% SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_id,null_id);

% Model Comparison Plots (Adjusted R-Squared)
SBJ04e_ERP_plot_RL_model_comparison_point(SBJ_id,an_id,stat_ids,null_id,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',0);
SBJ04e_ERP_plot_RL_model_comparison_R2_point(SBJ_id,an_id,stat_ids,null_id,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',1);
SBJ04e_ERP_plot_RL_model_comparison_R2_point(SBJ_id,an_id,stat_ids,null_id,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',0);
