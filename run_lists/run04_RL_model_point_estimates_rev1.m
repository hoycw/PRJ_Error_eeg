%% Reinforcement Learning based modeling and analyses of classic ERP metrics for Sequential PE Revision 1
% Developed over time, but last editted 4/16/21 by Colin W Hoy
% Final model_id = 'ERPEsL_DifFB', which includes all conditions
%       (no longer 'all', because that implies conditions are combined/averaged)
%   Sup. Fig. 3A: SBJ04e_ERP_plot_FRN_cond_metric_comparison_point
%   Sup. Fig. 3B and 3D: SBJ04c_ERP_grp_stats_LME_mean_window and SBJ04d_ERP_plot_stats_LME_mean_betas
%   Sup. Fig. 3C: SBJ04c_ERP_grp_stats_LME_P2P and SBJ04d_ERP_plot_stats_LME_p2p_betas
%   Sup. Fig. 3E/F/3G: SBJ04e_ERP_plot_RL_model_comparison_point
%   Sup. Fig. 9B: SBJ04c_ERP_p2p_latency_reg

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

%% ERP: Mean Window LME
an_id    = 'ERP_Fz_F2t1_dm2t0_fl05t20';% 'ERP_Pz_F2t1_dm2t0_fl05t20';

% Main RL Model in FRN window
stat_ids = {'VML_DifFB_lme_erpmn1FRN','SML_DifFB_lme_erpmn1FRN','ERPEsL_DifFB_lme_erpmn1FRN'};
% RL model in P3 window
% stat_ids = {'EsRPEL_DifFB_lme_erpmn1P3','ERPEs_DifFB_lme_erpmn1P3','ERPEsL_DifFB_lme_erpmn1P3'};%
% stat_ids = {'VML_DifFB_lme_erpmn1P3','SML_DifFB_lme_erpmn1P3','ERPEsL_DifFB_lme_erpmn1P3'};%
% RL model on Lik peak window
% stat_ids = {'EsRPEL_DifFB_lme_erpmn05Lik','ERPEs_DifFB_lme_erpmn05Lik','ERPEsL_DifFB_lme_erpmn05Lik'};%

plt_id    = 'bar_sigStar';
null_id   = 'SBJonly_all_lme_mn1FRN';

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

for st_ix = 1:numel(stat_ids)
    if ~isempty(strfind(stat_ids{st_ix},'erpmn'))
        % Average across ERPs (e.g., stat_id = 'ERPEsL_DifFB_lme_erpmn1FRN')
        SBJ04c_ERP_grp_stats_LME_mean_window(SBJ_id,proc_id,an_id,stat_ids{st_ix});
    else
        % Average across single trials (e.g., stat_id = 'ERPEsL_DifFB_lme_mn1FRN')
        %   Not used because literature typically averages over ERPs, not
        %   single trials, hence ERP mean
%         SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_id,stat_ids{st_ix});
    end
    
    % Sup. Fig. 3B and 3D: Plot mean window betas
    SBJ04d_ERP_plot_stats_LME_mean_betas(SBJ_id,proc_id,an_id,stat_ids{st_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Plot mean window betas with data as bar plot
%     SBJ04d_ERP_plot_stats_LME_mean_betas(SBJ_id,proc_id,an_id,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_data',1);
end

% Model Comparison Plots (AIC or Adjusted R-Squared)
%   Sup. Fig. 3E and 3G
SBJ04e_ERP_plot_RL_model_comparison_point(SBJ_id,an_id,stat_ids,null_id,plt_id,'AIC',save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype,'rm_null',0,'plot_null',0);
% SBJ04e_ERP_plot_RL_model_comparison_point(SBJ_id,an_id,stat_ids,null_id,plt_id,'R2',save_fig,...
%     'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',0);

%% ERP: Peak-to-Peak LME
% Main RL Model
an_id     = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_ids  = {'VML_DifFB_lme_p2pFRN','SML_DifFB_lme_p2pFRN','ERPEsL_DifFB_lme_p2pFRN'};% 
plot_peaks= 1;
plt_id    = 'bar_sigStar';

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for st_ix = 1:numel(stat_ids)
    % Compute peak-to-peak FRN metric and run RL model via LME
    SBJ04c_ERP_grp_stats_LME_P2P(SBJ_id,proc_id,an_id,stat_ids{st_ix},'plot_erps',1,'plot_peaks',1);
    
    % Sup. Fig. 3C: Plot RL model coefficients for peak-to-peak FRN metric
    SBJ04d_ERP_plot_stats_LME_p2p_betas(SBJ_id,an_id,stat_ids{st_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

    % Plot RL model coefficients and data from peak-to-peak FRN metric
%     SBJ04d_ERP_plot_stats_LME_p2p_betas(SBJ_id,an_id,stat_ids{st_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_data',1);%,'plot_latencies',1);
end

% Model Comparison Plots (Adjusted R-Squared) for FRN point estimates
%   Sup. Fig. 3F
SBJ04e_ERP_plot_RL_model_comparison_point(SBJ_id,an_id,stat_ids,'',plt_id,'AIC',save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',0);
% SBJ04e_ERP_plot_RL_model_comparison_point(SBJ_id,an_id,stat_ids,'',plt_id,'R2',save_fig,...
%     'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',0);

%% Peak Latency Regression
an_id     = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_ids  = {'ERPEsL_DifFB_lme_p2pFRN'};

proc_id   = 'eeg_full_ft';
save_fig  = 0;
fig_vis   = 'on';
fig_ftype = 'png';

for st_ix = 1:numel(stat_ids)
    % Run LME RL model on peak latencies from peak-to-peak FRN metric and plot latencies
%     SBJ_norm = 0;
%     SBJ04c_ERP_p2p_latency_reg(SBJ_id,proc_id,an_id,stat_ids{st_ix},SBJ_norm,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
%     % Sup. Fig. 9B: Run LME RL model on peak latencies from peak-to-peak FRN metric
%     %   also plot latencies after normalizing to mean peak latency within SBJ
    SBJ_norm = 1;
    SBJ04c_ERP_p2p_latency_reg(SBJ_id,proc_id,an_id,stat_ids{st_ix},SBJ_norm,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

stat_id  = 'ERPEsL_EHSu_lme_p2pFRN';
% Run LME RL model on peak latencies from peak-to-peak FRN metric and plot latencies
SBJ_norm = 0;
SBJ04c_ERP_p2p_latency_ttest(SBJ_id,an_id,stat_id,SBJ_norm);

%% FRN Metric Comparison: FRN by Condition
% Sup. Fig. 3A
an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
stat_ids  = {'ERPEsL_DifFB_lme_erpmn1FRN','ERPEsL_DifFB_lme_p2pFRN','ERPEsL_DifFB_lme_erpmn1P3'};

proc_id   = 'eeg_full_ft';
plt_id    = 'line_cond';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

% Sup. Fig. 1B: Plot mean window and peak-to-peak FRN point estimates
SBJ04e_ERP_plot_FRN_cond_metric_comparison_point(SBJ_id,proc_id,an_ids,stat_ids,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype);

% Plot peak-to-peak FRN and mean window (with mirrored y axis) point estimates for comparison
% SBJ04e_ERP_plot_FRN_cond_metric_comparison_point(SBJ_id,proc_id,an_id,stat_ids,plt_id,save_fig,...
%     'fig_vis',fig_vis,'fig_ftype',fig_ftype,'mirror_mean',1);


