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
proc_id  = 'eeg_full_ft';
stat_ids = {'ERPEsL_all_lme_st05t5'};

for s = 1:numel(SBJs)
    for st_ix = 1:numel(stat_ids)
        SBJ04a_RL_model(SBJs{s},proc_id,stat_ids{st_ix});
        % SBJ04b_BHV_RL_model_plot(SBJs{s},proc_id,stat_ids{st_ix});
    end
    close all;
end

%% ERP: Linear Mixed Effects Model (Over Time)
% Main RL Model
an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};%
stat_ids  = {'ERPEsL_all_lme_st05t5'};%'SML_all_lme_st05t5',
% stat_ids = {'S_all_lme_st05t5','V_all_lme_st05t5','sRPE_all_lme_st05t5'};
plt_id    = 'ts_F2t8_evnts_sigLine';
null_id   = 'SBJonly_all_lme_st05t5';

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
%       SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
      SBJ04d_ERP_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
%     SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_ids{an_ix},null_id);
    
    % Model Comparison Plots (Adjusted R-Squared)
%     SBJ04e_ERP_plot_RL_model_comparison_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',1);
%     SBJ04e_ERP_plot_RL_model_comparison_R2_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',1);
%     SBJ04e_ERP_plot_RL_model_comparison_R2_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',0);
end

%% ERP Topographies: Linear Mixed Effects Model (Mean Windows)
proc_id   = 'eeg_full_ft';
an_ids    = {'ERP_all_F2t1_dm2t0_fl05t20'};
stat_ids  = {'ERPEsL_all_lme_mn05sRPE','ERPEsL_all_lme_mn05uRPE','ERPEsL_all_lme_mn05Lik'};
% stat_ids  = {'ERPEsL_all_lme_mn05man1','ERPEsL_all_lme_mn05man225','ERPEsL_all_lme_mn05man325','ERPEsL_all_lme_mn05man375'};
% stat_ids  = {'VML_all_lme_mn1FRN','VML_all_lme_mn1P3'};
plt_id    = 'topo_F18t25';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
%         SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
%         SBJ04d_ERP_plot_stats_LME_RL_topo_reg(SBJ_id,an_ids{an_ix},stat_ids{st_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Plot Topo time series
    SBJ04d_ERP_plot_stats_LME_RL_topo_ts_reg(SBJ_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Power TFR: Linear Mixed Effects Model (Over Time-Frequency Power)
proc_id   = 'eeg_full_ft';
an_ids    = {'TFR_Fz_F2t1_db2t0_fl1t12','TFR_Pz_F2t1_db2t0_fl1t12'};%
stat_ids  = {'ERPEsL_all_lme_st0t5'};%'VML_all_lme_st0t5',
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

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
%         SBJ05d_PHS_grp_stats_wITPC_jkLME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
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

%% ========================================================================
%   OLD UNUSED ANALYSES (not going in the paper)
%  ========================================================================

% %% Pre-Feedback ERP: Linear Mixed Effects Model (Over Time)
% % Pre-Feedback RL Model
% an_ids    = {'ERP_Fz_F4t1_dm4t3_fl05t20'};%,'ERP_Pz_F4t1_dm4t3_fl05t20'};
% stat_ids = {'pWscr4D_all_lme_st3t5','rAscr4D_all_lme_st3t5'};
% % stat_ids = {'pWTaEr_all_lme_st3t0','pWThPr_all_lme_st3t0','pW4D_all_lme_st3t0','rATaEr_all_lme_st3t0','rAThPr_all_lme_st3t0','rA4D_all_lme_st3t0'};
% plt_id    = 'ts_F4t1_evnts_sigLine';%'ts_F2to1_evnts_sigLine';
% null_id   = 'SBJonly_all_lme_st3t5';%'SBJonly_all_lme_st3t5';
% 
% proc_id   = 'eeg_full_ft';
% save_fig  = 1;
% fig_vis   = 'on';
% fig_ftype = 'png';
% 
% % for an_ix = 1:numel(an_ids)
% %     for st_ix = 1:numel(stat_ids)
% %         SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
% %         SBJ04d_ERP_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
% %             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% %     end
% % %     SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_ids{an_ix},null_id);
% %     
% %     % Model Comparison Plots (Adjusted R-Squared)
% %     SBJ04e_ERP_plot_RL_model_comparison(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
% %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% % end
% 
% %% POW: Linear Mixed Effects Model (Over Time)
% % proc_id   = 'eeg_full_ft';
% % an_ids     = {'POW_Fz_F2t1_db2t0_fl4t8','POW_Fz_F2t1_db2t0_fl8t12','POW_Pz_F2t1_db2t0_fl1t4'};
% % %topo_plt_ids = {'topo_F18t25','topo_F18t25','topo_F3t45'};
% % stat_ids  = {'RLpRTulD_all_lme_st0t5'};%'RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5',
% % plt_id    = 'ts_F2to1_evnts_sigLine';
% % save_fig  = 1;
% % fig_vis   = 'on';
% % fig_ftype = 'png';
% % 
% % for an_ix = 1:numel(an_ids)
% %     for st_ix = 1:numel(stat_ids)
% %         SBJ05d_TFR_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
% %         SBJ05e_POW_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
% %             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% %     end
% %     
% %     % Model Comparison Plots (Adjusted R-Squared)
% % %     SBJ05f_POW_plot_RL_model_comparison(SBJ_id,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
% % %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% % end
