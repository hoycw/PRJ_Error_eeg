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
SBJ_id = 'goodOB';%'goodEEG1';
SBJs = fn_load_SBJ_list(SBJ_id);

%% Oddball ERP feature extraction
proc_id  = 'odd_full_ft';
feat_ids = {'N2sP3s_grpMW05'};
% Original options tested:
%--> OB ERP Features:
%       N2b = Novelty oddball N2
%       N2c = Target N2
%       P3a = Novelty oddball P3 (more central topo; max Cz)
%       P3b = Target P3 (more posterior topo; max Pz)
%--> Feature Metrics:
%       grpMW matches TT ERP window selection based on group
%       beta peaks and is also less noisy than individual peak
%       latencies/amplitudes that are vulnerable to overlapping components and SNR
%       1) grpMW- mean window around group peak latency
%       2) sbjMW- mean window around individual subject peak latency
%       3) sbjPk- peak amplitude of peak in individual subject ERP
%       4) p2p for N2: {'N2b_p2p','N2c_p2p'}
%           no results, likely suffers from low SNR (small N2) and component overlap confounds
%--> Window Length:
%       50 vs. 100 ms MW: using 50 ms to minimize component overlap
%--> Feature Combinations:
%       originally played with different groups of OB ERP features (e.g.,
%       N2bN2c, N2bP3a, etc.) to try to use multiple regression, but
%       they're all too correlated and N2sP3s is easier to run and allows
%       visualization of correlations in SBJ06b
%--> an_id is specified in the feat struct (always 'ERP_all_S2t1_dm2t0_fl05t20')

for ft_ix = 1:numel(feat_ids)
    % Extract Oddball ERP features
    if contains(feat_ids{ft_ix},'p2p')
        error('peak to peak is not reliable, so do not run this!');
%         SBJ06a_OB_ERP_save_p2p(SBJ_id,proc_id,feat_ids{ft_ix});
    else
%         SBJ06a_OB_ERP_save_mean_window(SBJ_id,proc_id,feat_ids{ft_ix});
    end
    
    % Plot correlations between ERP features before using as model predictors
    SBJ06b_OB_ERP_feature_corr(SBJ_id,proc_id,feat_ids{ft_ix});
end

%% Target Time ERP feature extraction
proc_id  = 'eeg_full_ft';
feat_ids = {'sRPE_Fz_grpMW05','uRPE_Pz_grpMW05','Lik_Fz_grpMW05'};
% Unused TT Features:
%   'Lik_Cz_grpMW05'
%   'P3_grpMW05','P3_sbjMW05'
%   'FRN_p2p','FRN_grpMW05','FRN_sbjMW05'
%   an_id in ft struct {'ERP_Fz_F2t1_dm2t0_fl05t20', 'ERP_Pz_F2t1_dm2t0_fl05t20'}

for ft_ix = 1:numel(feat_ids)
    % Extract Target Time ERP features
    if contains(feat_ids{ft_ix},'p2p')
        error('peak to peak is not reliable, so do not run this!');
%         SBJ06c_TT_ERP_save_p2p(SBJ_id,proc_id,feat_ids{ft_ix});
    else
        SBJ06c_TT_ERP_save_mean_window(SBJ_id,proc_id,feat_ids{ft_ix});
    end
end

%% Oddball vs. Target Time ERP Comparison: Point Estimates
tt_proc_id = 'eeg_full_ft';
ob_proc_id = 'odd_full_ft';
model_id   = 'ERPEsL_DifFB';

% FRN Parameters:
stat_ids   = {...
    'N2sP3s_grpMW05_DifFB_corr_sRPE_Fz_grpMW05'...
    'N2sP3s_grpMW05_DifFB_corr_uRPE_Pz_grpMW05'...
    'N2sP3s_grpMW05_DifFB_corr_Lik_Fz_grpMW05'...
    };

% plt_id    = 'bar_sigStar';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for st_ix = 1:numel(stat_ids)
%     SBJ06d_OB_TT_ERP_grp_stats_corr_pt(SBJ_id,tt_proc_id,ob_proc_id,stat_ids{st_ix},...
%         'save_fig',save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    SBJ06e_OB_TT_ERP_grp_plot_corr_reg_comparison(SBJ_id,tt_proc_id,stat_ids{st_ix},model_id,...
        'save_fig',save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Oddball TFR feature extraction
proc_id  = 'odd_full_ft';
feat_ids = {'deltaRare_2t5','thetaRare_2t4'};%,'deltaRare_3t5'};
    %,'thetaOdd_2t4','thetaTar_2t4',...'deltaOdd_3t5','deltaTar_3t5'};
%   an_id is specified in the feat struct (always 'ERP_all_S2t1_dm2t0_fl05t20')

for ft_ix = 1:numel(feat_ids)
    % Extract Oddball TFR features
    SBJ08a_OB_TFR_save_mean_window(SBJ_id,proc_id,feat_ids{ft_ix});
    
    % If multiple TFR features, plot correlations before using as model predictors
    if contains(feat_ids{ft_ix},'Rare')
        SBJ08b_OB_TFR_feature_corr(SBJ_id,proc_id,feat_ids{ft_ix});
    end
end

%% Target Time TFR feature extraction
proc_id  = 'eeg_full_ft';
feat_ids = {'thetaFRN_2t4','deltaP3_2t45'};%'thetaFRN_2t35','deltaP3_3t45'};
%   an_id is specified in the feat struct (always 'ERP_all_S2t1_dm2t0_fl05t20')

for ft_ix = 1:numel(feat_ids)
    % Extract Oddball TFR features
    SBJ08c_TT_TFR_save_mean_window(SBJ_id,proc_id,feat_ids{ft_ix});
end

%% OB-TT TFR comparison
tt_proc_id = 'eeg_full_ft';
ob_proc_id = 'odd_full_ft';
stat_ids   = {'deltaRare_2t5_DifFB_corr_deltaP3_2t45'};%'thetaRare_2t4_DifFB_corr_thetaFRN_2t4'};
%'thetaRare_2t4_DifFB_corr_thetaFRN_2t35','deltaRare_3t5_DifFB_corr_deltaP3_2t45'};
    
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';
for st_ix = 1:numel(stat_ids)
    SBJ08d_OB_TT_TFR_grp_stats_corr_pt(SBJ_id,tt_proc_id,ob_proc_id,...
        stat_ids{st_ix},'save_fig',save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end
