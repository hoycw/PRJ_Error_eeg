%% Reinforcement Learning based modeling and analyses for Sequential PE Initial Submission
% Developed over time, but editted 8/X/20 by Colin W Hoy
% Final model_id = 'ERPEsL_all'
%   Fig. 1D: SBJ04b_BHV_RL_model_plot and SBJ04b_BHV_RL_model_plot_grp
%   Fig. 2: 

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
SBJ_id = 'goodall';%'good1';%'good2';%
SBJs = fn_load_SBJ_list(SBJ_id);

%% Single SBJ RL Model
proc_id   = 'eeg_full_ft';
stat_ids  = {'ERPEsL_all_lme_st05t5'};
% Alternative (worse) models: 'RSVPE_all_lme_mn1FRN','SML_all_lme_mn1FRN','VML_all_lme_mn1FRN'
fig_vis   = 'on';
save_fig  = 1;
fig_ftype = 'png';

for s = 1:numel(SBJs)
    for st_ix = 1:numel(stat_ids)
        % Run model
        SBJ04a_RL_model(SBJs{s},proc_id,stat_ids{st_ix});
        
        % Fig. 1D: Plot model fit to tolerance and outcomes/accuracy
        SBJ04b_BHV_RL_model_plot(SBJs{s},proc_id,stat_ids{st_ix},...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    close all;
end

% Fig. 1D inset: Plot group model fits (overlapping sigmoids without tolerance)
SBJ04b_BHV_RL_model_plot_grp(SBJ_id,proc_id,stat_id,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype);

% Sup. Fig. 1A: Plot model predicitons by condition across group
plt_id    = 'line_cond';
for st_ix = 1:numel(stat_ids)
    SBJ04b_plot_model_predictions(SBJ_id,proc_id,stat_ids{st_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% ERP: Linear Mixed Effects Model (Over Time)
% Main RL Model
an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
stat_ids  = {'ERPEsL_all_lme_st05t5'};
plt_id    = 'ts_F2t8_evnts_sigLine';
null_id   = 'SBJonly_all_lme_st05t5';

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
      SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
      SBJ04d_ERP_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    % Optional: run SBJ only model (random intercepts, no regressors) for
    % baseline model performance; not in the paper
%     SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_ids{an_ix},null_id);

    % Model Comparison Plots (Adjusted R-Squared)
%     SBJ04e_ERP_plot_RL_model_comparison_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',1);
%     SBJ04e_ERP_plot_RL_model_comparison_R2_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',1);
%     SBJ04e_ERP_plot_RL_model_comparison_R2_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',0);
end

% Electrode R2 Comparison Plot
plt_id     = 'ts_F0t5_evnts_sigLine';
for st_ix = 1:numel(stat_ids)
%     SBJ04e_ERP_plot_RL_elec_comparison_R2_ts(SBJ_id,an_ids,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted');
end

%% Beta Topographies: Linear Mixed Effects Model (Mean Windows)
proc_id   = 'eeg_full_ft';
an_ids    = {'ERP_all_F2t1_dm2t0_fl05t20'};
stat_ids  = {'ERPEsL_all_lme_mn05sRPE','ERPEsL_all_lme_mn05uRPE','ERPEsL_all_lme_mn05Lik'};
% stat_ids  = {'ERPEsL_all_lme_mn05man1','ERPEsL_all_lme_mn05man225','ERPEsL_all_lme_mn05man325','ERPEsL_all_lme_mn05man375'};'ERPEsL_all_lme_mn150man375'};%
% stat_ids  = {'VML_all_lme_mn1FRN','VML_all_lme_mn1P3'};
plt_id    = 'topo_F18t25';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ04d_ERP_plot_stats_LME_RL_topo_reg(SBJ_id,an_ids{an_ix},stat_ids{st_ix},...
            plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Plot Topo time series
%     SBJ04d_ERP_plot_stats_LME_RL_topo_ts_reg(SBJ_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);

    % ERP Topographies: ERP topo dynamics in model mean windows
    conditions = 'DifFB';
    SBJ03c_ERP_plot_grp_topo_ts_cond(SBJ_id,conditions,proc_id,an_ids{an_ix},stat_ids,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Power TFR: Linear Mixed Effects Model (Over Time-Frequency Power)
proc_id   = 'eeg_full_ft';
an_ids    = {'TFR_Fz_F2t1_db2t0_fl1t12','TFR_Pz_F2t1_db2t0_fl1t12'};%
stat_ids  = {'ERPEsL_all_lme_st0t5'};%'VML_all_lme_st0t5',
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
%         SBJ05d_TFR_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ05e_TFR_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
end

%% PHS TFR: Linear Mixed Effects Model (Over Time-Frequency Phase)
% conditions = 'DifFB';
proc_id   = 'eeg_full_ft';
an_ids    = {'PHS_Fz_F2t1_fl1t12'};%,'PHS_Pz_F2t1_fl1t12'
model_ids = {'ERPEsL_all'};
model_win = 'st0t5';
%plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(model_ids)
%         % Circular-Linear Regression: Group level, separate for each regressor
        stat_id = [model_ids{st_ix} '_CLreg_' model_win];
%         SBJ05d_PHS_grp_stats_CLreg_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
        SBJ05e_PHS_plot_stats_CLreg_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Jack-Knife wITPC LME Regression: Group level
        stat_id = [model_ids{st_ix} '_lme_' model_win];
%         SBJ05d_PHS_grp_stats_ITPC_jkLME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
        SBJ05e_TFR_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        
%         % wITPC: SBJ level z-score, GRP t-test
%         stat_id = [model_ids{st_ix} '_wITPC_' model_win];
% %         SBJ05d_PHS_grp_stats_wITPC_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
% %         SBJ05e_PHS_plot_stats_wITPC_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
% %             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% 
% %         % Circular-Linear Correlation: Group Level
%         stat_id = [model_ids{st_ix} '_CLcorr_' model_win];
% %         SBJ05d_PHS_grp_stats_CLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
% %         SBJ05e_PHS_plot_stats_CLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
% %             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         
% %         % Circular-Linear Correlation: SBJ z-score, GRP t-test
% %         stat_id = [model_ids{st_ix} '_zCLcorr_' model_win];
% %         SBJ05d_PHS_grp_stats_zCLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
% %         SBJ05e_PHS_plot_stats_zCLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
% %             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
    end
    
    % Model Comparison Plots (Adjusted R-Squared)
%     error('this needs to be a martix version!');
%     SBJ05f_TFR_plot_RL_model_comparison(SBJ_id,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end


