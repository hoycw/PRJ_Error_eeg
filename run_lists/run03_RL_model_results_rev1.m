%% Reinforcement Learning based modeling and analyses for Sequential PE Initial Submission
% Developed over time, but editted 8/X/20 by Colin W Hoy
% Final model_id = 'ERPEsL_all' (all = 'DifFB', which includes all conditions)
%   Fig. 1D: SBJ04b_BHV_RL_model_plot and SBJ04b_BHV_RL_model_plot_grp
%   Fig. 2: SBJ04c_ERP_grp_stats_LME_RL and SBJ04d_ERP_plot_stats_LME_RL_fits
%   Fig. 3: SBJ04c_ERP_grp_stats_LME_RL and SBJ04d_ERP_plot_stats_LME_RL_topo_ts_reg
%   Fig. 4: SBJ05d_TFR_grp_stats_LME_RL and SBJ05e_TFR_plot_stats_LME_RL_fits
%   Sup. Fig. 1A: SBJ04b_plot_model_predictions
%   Sup. Fig. 2: SBJ04e_ERP_plot_RL_elec_comparison_R2_ts
%   Sup. Fig. 3: SBJ03c_ERP_plot_grp_topo_ts_cond
%   Sup. Fig. 6C: SBJ05d_PHS_grp_stats_CLreg_RL and SBJ05e_PHS_plot_stats_CLreg_RL
%   Sup. Fig. 6D: SBJ05d_PHS_grp_stats_ITPC_jkLME_RL and SBJ05e_TFR_plot_stats_LME_RL_fits

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
SBJ_id = 'goodall';%'good2';%'good1';%
SBJs = fn_load_SBJ_list(SBJ_id);

%% Single SBJ RL Model
proc_id   = 'eeg_full_ft';
stat_ids  = {'EuRPEL_Pos_lme_st05t5'};%'ERPEsL_DifFB_lme_st05t5'};

% Subjective rating bias models:
% stat_ids  = {'ERBuRPE_EzSu_lme_st05t5'};%'ERPEsL_pW25hd_DifFB_lme_st05t5','ERPEsL_pW25_DifFB_lme_st05t5'};
% Auditory Salience models:
% stat_ids  = {'rsRPE_EHSu_lme_st05t5'};%'ERBsRPE_EHSu_lme_st05t5','ERBrsRPE_EHSu_lme_st05t5'};%'ERB_EHSu_lme_st0t5','ERBr_EHSu_lme_st0t5','rough_EHSu_lme_st0t5','AudSal_EHSu_lme_st0t5'};
% Outcome-based models:
% stat_ids  = {'VML_DifFB_lme_st05t5','SML_DifFB_lme_st05t5'};%'ERPEs_DifFB_lme_st05t5'};%'uRPE_Neg_lme_st05t5'};%'ML_Neg_lme_st05t5'};%'ERPEsL_all_lme_st05t5'};
%   for plotting valence and value predictions: 'VSML_DifFB_lme_st05t5'

fig_vis   = 'on';
save_fig  = 1;
fig_ftype = 'png';

for s = 2:numel(SBJs)
    for st_ix = 1:numel(stat_ids)
        % Run model
        SBJ04a_RL_model(SBJs{s},proc_id,stat_ids{st_ix});
        
        % Fig. 1D: Plot model fit to tolerance and outcomes/accuracy
%         SBJ04b_BHV_RL_model_plot(SBJs{s},proc_id,stat_ids{st_ix},...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    close all;
end

% Fig. 1D inset: Plot group model fits (overlapping sigmoids without tolerance)
% SBJ04b_BHV_RL_model_plot_grp(SBJ_id,proc_id,stat_id,...
%     'fig_vis',fig_vis,'fig_ftype',fig_ftype);

% Sup. Fig. 1A: Plot model predicitons by condition across group
plt_id    = 'line_cond';
for st_ix = 1:numel(stat_ids)
%     SBJ04b_RL_model_plot_grp_predictions(SBJ_id,proc_id,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% ERP: Linear Mixed Effects Model (Over Time)
% Plots Fig. 2 and 3; Sup. Fig. 2
% Main RL Model
an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};%,'ERP_Cz_F2t1_dm2t0_fl05t20'};%,

% Main RL model:
% stat_ids  = {'ERPEsL_DifFB_lme_st05t5'};
% Split negative and positive outcomes:
% stat_ids  = {'uRPEL_Neg_lme_st05t5','uRPEL_Pos_lme_st05t5'};%'uRPE_Neg_lme_st05t5','uRPE_Pos_lme_st05t5'};
% hard subjective bias:
stat_ids  = {'ERPEsL_pW25hd_DifFB_lme_st05t5'};
% Auditory Salience:
% stat_ids  = {'ERBsRPE_EHSu_lme_st05t5'};
% Outcome models:
% stat_ids  = {'VML_DifFB_lme_st05t5','SML_DifFB_lme_st05t5'};%,'ERPEsL_DifFB_lme_st05t5'};
% RL models without uRPE or Lik:
% stat_ids  = {'EsRPEL_DifFB_lme_st05t5','ERPEs_DifFB_lme_st05t5'};%,'ERPEsL_DifFB_lme_st05t5'};
% All models:
% stat_ids  = {'VML_DifFB_lme_st05t5','SML_DifFB_lme_st05t5','EsRPEL_DifFB_lme_st05t5','ERPEs_DifFB_lme_st05t5','ERPEsL_DifFB_lme_st05t5'};

plt_id    = 'ts_F2t8_evnts_sigLine';%'ts_F0t6_evnts_sigLine';
null_id   = 'SBJonly_all_lme_st05t5';

proc_id   = 'eeg_full_ft';
save_fig  = 0;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
      % Run LME RL model on ERPs over time
      SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
      
%       % Fig. 2: Plot model results (ERPs, coefficients, model fit)
%       SBJ04d_ERP_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);

      % Fig. 2ABCD: Plot model results (ERPs, coefficients, model fit)
      %     with window overlays for N3, P3, sRPE, uRPE, Lik
      if contains(an_ids{an_ix},'Fz')
          erp_win_name = {'N2'};
          if contains(stat_ids{st_ix},'ERPEsL');  beta_win_name = {'sRPE','Lik'};
          elseif contains(stat_ids{st_ix},'VML'); beta_win_name = {'Val','Lik'};
          elseif contains(stat_ids{st_ix},'SML'); beta_win_name = {'Sign','Lik'}; end
      elseif contains(an_ids{an_ix},'Pz')
          erp_win_name = {'P3'};
          if contains(stat_ids{st_ix},'EsRPE'); beta_win_name = {'Lik'};
          else                                  beta_win_name = {'uRPE'}; end
      end
%       SBJ04d_ERP_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype,'erp_win_name',erp_win_name,'beta_win_name',beta_win_name);
    end
    
    % Optional: run SBJ only model (random intercepts, no regressors) for
    %   baseline model performance; not in the paper
%     SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_ids{an_ix},null_id);

    % Model Comparison Plots: AIC Performance Relative to SBJonly null model
    if contains(an_ids{an_ix},'Fz')
%         aic_mean_reg = {'sRPE','ERPEsL_DifFB_lme_st05t5'};
        aic_mean_win = [-0.025 0.025]+0.216;    % Peak beta for sRPE in goodall
    elseif contains(an_ids{an_ix},'Pz')
% % %         aic_mean_reg = {'uRPE','ERPEsL_DifFB_lme_st05t5'};
        aic_mean_win = [-0.025 0.025]+0.308;    % beta peaks in goodall for uRPE (0.308) and Lik (0.38)
    end
%     aic_mean_win = [-0.025 0.025]+0.380;
%     SBJ04e_ERP_plot_RL_model_comparison_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'rm_null',1,'mean_win',aic_mean_win);%'mean_reg',aic_mean_reg);

    % Model Comparison Plots: AIC Performance
%     SBJ04e_ERP_plot_RL_model_comparison_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',1);%,'mean_reg',aic_mean_reg);

    % Print model comparisons: AIC Performance
    aic_mean_win = [-0.025 0.025]+[0.216 0.308 0.380]';
    for mw_ix = 1:size(aic_mean_win,1)
        SBJ04e_ERP_print_RL_model_comparison_win(SBJ_id,an_ids{an_ix},stat_ids,null_id,...
            aic_mean_win(mw_ix,:),'rm_null',1);
    end

    % Model Comparison Plots: R2 Fits Relative to SBJonly null model
%     SBJ04e_ERP_plot_RL_model_comparison_R2_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',1);

    % Model Comparison Plots: R2 Fits Overall
%     SBJ04e_ERP_plot_RL_model_comparison_R2_ts(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',0);
end

% Sup. Fig. 2: Electrode R2 Comparison Plot
plt_id     = 'ts_F0t5_evnts_sigLine';   % plot only for stat_lim
for st_ix = 1:numel(stat_ids)
%     SBJ04e_ERP_plot_RL_elec_comparison_R2_ts(SBJ_id,an_ids,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted');
end

%% Beta Topographies: Linear Mixed Effects Model (Mean Windows)
% Plots Fig. 3 and Sup. Fig. 3
proc_id   = 'eeg_full_ft';
an_ids    = {'ERP_all_F2t1_dm2t0_fl05t20'};
stat_ids  = {'uRPEL_Neg_lme_mn05man216','uRPEL_Neg_lme_mn05man308','uRPEL_Neg_lme_mn05man380'};
% stat_ids  = {'uRPEL_Pos_lme_mn05man216','uRPEL_Pos_lme_mn05man308','uRPEL_Pos_lme_mn05man380'};
% stat_ids  = {'uRPE_Neg_lme_mn05man216','uRPE_Pos_lme_mn05man216','uRPE_Neg_lme_mn05man308','uRPE_Pos_lme_mn05man308'};
% stat_ids = {'ERPEsL_all_lme_mn05sRPE','ERPEsL_all_lme_mn05uRPE','ERPEsL_all_lme_mn05Lik'};
plt_id    = 'topo_F18t25';
save_fig  = 0;
fig_vis   = 'on';
fig_ftype = 'svg';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        % Run LME RL model on ERPs averaged in time window for all electrodes
%         SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
        
        % Plot individual model coefficient topographies
%         SBJ04d_ERP_plot_stats_LME_RL_topo_reg(SBJ_id,an_ids{an_ix},stat_ids{st_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Fig. 3: Plot Beta Topo across time points
    SBJ04d_ERP_plot_stats_LME_RL_topo_ts_reg(SBJ_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Sup. Fig. 3: ERP topography dynamics in model averaging windows
    conditions = 'DifFB';
%     SBJ03c_ERP_plot_grp_topo_ts_cond(SBJ_id,conditions,proc_id,an_ids{an_ix},stat_ids,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Power TFR: Linear Mixed Effects Model (Over Time-Frequency Power)
% Plots Fig. 4 and Fig. 6
an_ids    = {'TFR_Fz_F2t1_db2t0_fl1t12','TFR_Pz_F2t1_db2t0_fl1t12'};
% Main RL Model (Fig. 4):
% stat_ids  = {'ERPEsL_all_lme_st0t5'};

% Split Positive and Negative Conditions (Fig. 6):
stat_ids  = {'uRPEL_Neg_lme_st0t5','uRPEL_Pos_lme_st0t5'};

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        % Run LME RL model on TFR power
%         SBJ05d_TFR_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
        
        % Fig. 4: Plot TFR Power model coefficients
        SBJ05e_TFR_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
end

%% PHS TFR: Linear Mixed Effects Model (Over Time-Frequency Phase)
% Plots Sup. Fig. 6C and 6D
proc_id   = 'eeg_full_ft';
an_ids    = {'PHS_Fz_F2t1_fl1t12'};%,'PHS_Pz_F2t1_fl1t12'
model_ids = {'uRPEL_Neg','uRPEL_Pos'};%'ERPEsL_all'};
model_win = 'st0t5';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(model_ids)
        % Circular-Linear Regression: Group level, separate for each regressor
        stat_id = [model_ids{st_ix} '_CLreg_' model_win];
        SBJ05d_PHS_grp_stats_CLreg_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
        
        % Sup. Fig. 6C: Plot C-L Regression TFR Phase Coefficients
        SBJ05e_PHS_plot_stats_CLreg_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Jack-Knife ITPC LME Regression: Group level
        stat_id = [model_ids{st_ix} '_lme_' model_win];
        SBJ05d_PHS_grp_stats_ITPC_jkLME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
        
        % Sup. Fig. 6D: Plot ITPC TFR Phase Coefficients
        SBJ05e_TFR_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
end


