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

