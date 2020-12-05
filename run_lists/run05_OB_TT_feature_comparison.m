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
feat_ids = {'N2bN2c_grpMW05','N2bN2c_sbjMW05'};%'N2bP3a_sbjMW05','N2bP3a_grpMW05','N2cP3b_sbjMW05','N2cP3b_grpMW05'};
    %'P3aP3b_grpMW1'};%'N2cP3b_sbjPk','N2cP3b_sbjMW1','N2cP3b_grpMW1'};
    %'N2sP3s_grpMW05','N2sP3s_sbjMW05','N2sP3s_grpMW1','N2sP3s_sbjMW1','N2sP3s_sbjPk'};
% feat_ids = {'N2bP3a_sbjPk','N2bP3a_sbjMW1','N2bP3a_grpMW1'};
%   an_id is specified in the feat struct (always 'ERP_all_S2t1_dm2t0_fl05t20')

for ft_ix = 1:numel(feat_ids)
    % Extract Oddball ERP features
    if contains(feat_ids{ft_ix},'p2p')
        error('not ready yet');
    else
        SBJ06a_OB_ERP_save_mean_window(SBJ_id,proc_id,feat_ids{ft_ix});
    end
    
    % Plot correlations between ERP features before using as model predictors
    SBJ06b_OB_ERP_feature_corr(SBJ_id,proc_id,feat_ids{ft_ix});
end

%% Target Time ERP feature extraction
proc_id  = 'eeg_full_ft';
feat_ids = {'FRN_grpMW05','FRN_sbjMW05'};
%   an_id in ft struct {'ERP_Fz_F2t1_dm2t0_fl05t20', 'ERP_Pz_F2t1_dm2t0_fl05t20'}

for ft_ix = 1:numel(feat_ids)
    % Extract Target Time ERP features
    if contains(feat_ids{ft_ix},'p2p')
        error('not ready yet');
    else
        SBJ06c_TT_ERP_save_mean_window(SBJ_id,proc_id,feat_ids{ft_ix});
    end
end

%% Oddball vs. Target Time ERP Comparison
tt_proc_id = 'eeg_full_ft';
ob_proc_id = 'odd_full_ft';

% FRN Parameters:
stat_ids   = {...
    'N2bP3a_grpMW05_DifFB_reg_FRNgrpMW05','N2bP3a_sbjMW05_DifFB_reg_FRNsbjMW05'};
%     'N2cP3b_grpMW05_DifFB_reg_FRNgrpMW05','N2cP3b_sbjMW05_DifFB_reg_FRNsbjMW05'};

% Original stuff:
%     'P3aP3b_grpMW1_All_reg_mn05LikCz'};
%     'N2cP3b_grpMW1_AllNeg_reg_erpmn1FRN', 'N2cP3b_sbjMW1_AllNeg_reg_erpmn1FRN'...
%     'N2cP3b_grpMW1_AllPos_reg_erpmn1FRN', 'N2cP3b_sbjMW1_AllPos_reg_erpmn1FRN'...
%     };
% stat_ids   = {...
%     'N2bP3a_grpMW1_All_reg_erpmn1FRN', 'N2bP3a_sbjMW1_All_reg_erpmn1FRN',...
%     'N2cP3b_grpMW1_All_reg_erpmn1FRN', 'N2cP3b_sbjMW1_All_reg_erpmn1FRN'...
%     };
% P3 Parameters:
% an_id      = 'ERP_Pz_F2t1_dm2t0_fl05t20';
% stat_ids   = {...
%     'P3aP3b_grpMW1_All_reg_mn05uRPE'...
%     'N2bP3a_grpMW1_All_reg_erpmn1P3', 'N2bP3a_sbjMW1_All_reg_erpmn1P3'...
%     'N2cP3b_grpMW1_All_reg_erpmn1P3', 'N2cP3b_sbjMW1_All_reg_erpmn1P3',...
%     };

% plt_id    = 'bar_sigStar';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for st_ix = 1:numel(stat_ids)
    SBJ06d_OB_TT_ERP_grp_stats_reg(SBJ_id,tt_proc_id,ob_proc_id,...
        stat_ids{st_ix},'save_fig',save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

