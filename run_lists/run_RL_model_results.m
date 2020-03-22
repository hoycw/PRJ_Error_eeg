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
sbj_file = fopen([root_dir 'PRJ_Error_EEG/scripts/SBJ_lists/' SBJ_id '.sbj']);
tmp = textscan(sbj_file,'%s');
fclose(sbj_file);
SBJs = tmp{1}; clear tmp;

%% Single SBJ RL Model
proc_id  = 'eeg_full_ft';
% stat_ids = {'RLbA_all_lme_st0t5','RLbApW_all_lme_st0t5','RLrA_all_lme_st0t5','RLrApW_all_lme_st0t5'};
stat_ids  = {'RLTaEr_all_lme_st1t5','RLThPr_all_lme_st1t5','RLsD_all_lme_st1t5','RLuD_all_lme_st1t5','RLsTauTh_all_lme_st1t5','RLuTasTh_all_lme_st1t5'};
% stat_ids  = {'RL_all_lme_st0t5','RL3D_all_lme_st0t5','RL4D_all_lme_st0t5'};
% stat_ids = {'rATar_all_lme_st3t5','rAallD_all_lme_st3t5'};
% stat_ids = {'pWTaEr_all_lme_st3t5','pWThPr_all_lme_st3t5','pW4D_all_lme_st3t5','rATaEr_all_lme_st3t5','rAThPr_all_lme_st3t5','rA4D_all_lme_st3t5'};
% RL models:
%   RL3D: pWin, sPE, uPE, sTar, uTar, uThr (no sThr!!!)
%   RLfullD: pWin, sPE, uPE, sTar, uTar, uThr, sThr
%   RL: pWin, sPE, uPE
%   pWTarD: pWin, uTar, sTar
%   pWallD: pWin, uTar, sTar, uThr, sThr

for s = 1:numel(SBJs)
    for st_ix = 1:numel(stat_ids)
        SBJ04a_RL_model(SBJs{s},proc_id,stat_ids{st_ix});
        % SBJ04b_BHV_RL_model_plot(SBJs{s},proc_id,stat_ids{st_ix});
    end
    close all;
end

%% ERP: Linear Mixed Effects Model (Over Time)
% Main RL Model
%an_ids    = {'ERPlp_Fz_F2t1_dm2t0_fl05t20'};
an_ids    = {'ERP_Pz_F2t1_dm2t0_fl05t20'};%,'ERP_Pz_F2t1_dm2t0_fl05t20'};
stat_ids  = {'RL3D_all_lme_st1t5'};%,'RLTaEr_all_lme_st1t5','RLThPr_all_lme_st1t5','RLsD_all_lme_st1t5','RLuD_all_lme_st1t5','RLsTauTh_all_lme_st1t5','RLuTasTh_all_lme_st1t5'};
% stat_ids  = {'RL_all_lme_st0t5','RL3D_all_lme_st0t5','RL4D_all_lme_st0t5','RLallD_all_lme_st0t5'};
% stat_ids = {'RLbA_all_lme_st0t5','RLrA_all_lme_st0t5','RL_all_lme_st0t5','RLrApW_all_lme_st0t5'};%'RLrApW_all_lme_st0t5'};
% stat_ids  = {'RL_all_lme_st0t5','RL3D_all_lme_st0t5','RLfullD_all_lme_st0t5'};
plt_id    = 'ts_F2to1_evnts_sigLine';
null_id   = 'SBJonly_all_lme_st1t5';

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
%     SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_ids{an_ix},null_id);
    
    % Model Comparison Plots (Adjusted R-Squared)
%     SBJ04e_ERP_plot_RL_model_comparison(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_null',1);
%     SBJ04e_ERP_plot_RL_model_comparison_R2(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',1);
%     SBJ04e_ERP_plot_RL_model_comparison_R2(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'r2_version','Adjusted','rm_null',0);
end

%% Pre-Feedback ERP: Linear Mixed Effects Model (Over Time)
% Pre-Feedback RL Model
an_ids    = {'ERP_Fz_F4t1_dm4t3_fl05t20'};%,'ERP_Pz_F4t1_dm4t3_fl05t20'};
stat_ids = {'pWscr4D_all_lme_st3t5','rAscr4D_all_lme_st3t5'};
% stat_ids = {'pWTaEr_all_lme_st3t0','pWThPr_all_lme_st3t0','pW4D_all_lme_st3t0','rATaEr_all_lme_st3t0','rAThPr_all_lme_st3t0','rA4D_all_lme_st3t0'};
plt_id    = 'ts_F4t1_evnts_sigLine';%'ts_F2to1_evnts_sigLine';
null_id   = 'SBJonly_all_lme_st3t5';%'SBJonly_all_lme_st3t5';

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
%     SBJ04c_ERP_grp_stats_LME_SBJonly(SBJ_id,proc_id,an_ids{an_ix},null_id);
    
    % Model Comparison Plots (Adjusted R-Squared)
    SBJ04e_ERP_plot_RL_model_comparison(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% ERP: Linear Mixed Effects Model (Mean Windows)
proc_id   = 'eeg_full_ft';
an_ids    = {'ERP_all_F2t1_dm2t0_fl05t20'};%'ERPlp_all_F2t1_dm2t0_fl05t20'};
stat_ids  = {'RL3D_all_lme_mn05sPE05','RL3D_all_lme_mn05uPE05'};
% stat_ids  = {'RLpRTlD_all_lme_mn2t3','RLpRTlD_all_lme_mn3t4'};
plt_ids   = {'topo_F18t25','topo_F3t45'};
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        SBJ04c_ERP_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ04d_ERP_plot_stats_LME_RL_topo_reg(SBJ_id,an_ids{an_ix},stat_ids{st_ix},...
            plt_ids{st_ix},save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Model Comparison Plots (Adjusted R-Squared)
%     error('need topo version!');
%     SBJ04e_ERP_plot_RL_model_comparison(SBJ_id,an_ids{an_ix},stat_ids,null_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Power TFR: Linear Mixed Effects Model (Over Time-Frequency Power)
% conditions = 'DifFB';
proc_id   = 'eeg_full_ft';
an_ids     = {'TFR_Fz_F2t1_db2t0_fl1t12b05','TFR_Pz_F2t1_db2t0_fl1t12b05'};
stat_ids  = {'RL3D_all_lme_st0t5'};
%plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        SBJ05d_TFR_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ05e_TFR_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Model Comparison Plots (Adjusted R-Squared)
%     error('this needs to be a martix version!');
%     SBJ05f_TFR_plot_RL_model_comparison(SBJ_id,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% PHS TFR: Linear Mixed Effects Model (Over Time-Frequency Phase)
% conditions = 'DifFB';
proc_id   = 'eeg_full_ft';
an_ids    = {'PHS_Fz_F2t1_fl1t12b05'};%,'TFR_Pz_F2t1_db2t0_fl1t12b05'};
model_ids = {'RL_all'};
model_win = 'st0t5';
%stat_ids  = {'RLpRTulD_all_lme_st0t5'};%'RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5',
%plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(model_ids)
        % Circular-Linear Correlation: Group Level
        stat_id = [model_ids{st_ix} '_CLcorr_' model_win];
%         SBJ05d_PHS_grp_stats_CLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
        SBJ05e_PHS_plot_stats_CLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
%         % Circular-Linear Correlation: SBJ z-score, GRP t-test
%         stat_id = [model_ids{st_ix} '_zCLcorr_' model_win];
%         SBJ05d_PHS_grp_stats_zCLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
%         SBJ05e_PHS_plot_stats_zCLcorr_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
%         % Circular-Linear Regression: Group level, separate for each regressor
%         stat_id = [model_ids{st_ix} '_CLreg_' model_win];
%         SBJ05d_PHS_grp_stats_CLreg_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
%         SBJ05e_PHS_plot_stats_CLreg_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
%         % Jack-Knife wITPC LME Regression: Group level
%         stat_id = [model_ids{st_ix} '_lme_' model_win];
%         SBJ05d_PHS_grp_stats_wITPC_jkLME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
%         SBJ05e_TFR_plot_stats_LME_RL_fits(tfr_SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
%         % wITPC: SBJ level z-score, GRP t-test
%         stat_id = [model_ids{st_ix} '_wITPC_' model_win];
%         SBJ05d_PHS_grp_stats_wITPC_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id);
%         SBJ05e_PHS_plot_stats_wITPC_RL(SBJ_id,proc_id,an_ids{an_ix},stat_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Model Comparison Plots (Adjusted R-Squared)
%     error('this needs to be a martix version!');
%     SBJ05f_TFR_plot_RL_model_comparison(SBJ_id,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% POW: Linear Mixed Effects Model (Over Time)
% proc_id   = 'eeg_full_ft';
% an_ids     = {'POW_Fz_F2t1_db2t0_fl4t8','POW_Fz_F2t1_db2t0_fl8t12','POW_Pz_F2t1_db2t0_fl1t4'};
% %topo_plt_ids = {'topo_F18t25','topo_F18t25','topo_F3t45'};
% stat_ids  = {'RLpRTulD_all_lme_st0t5'};%'RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5',
% plt_id    = 'ts_F2to1_evnts_sigLine';
% save_fig  = 1;
% fig_vis   = 'on';
% fig_ftype = 'png';
% 
% for an_ix = 1:numel(an_ids)
%     for st_ix = 1:numel(stat_ids)
%         SBJ05d_TFR_grp_stats_LME_RL(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix});
%         SBJ05e_POW_plot_stats_LME_RL_fits(SBJ_id,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     end
%     
%     % Model Comparison Plots (Adjusted R-Squared)
% %     SBJ05f_POW_plot_RL_model_comparison(SBJ_id,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
% %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% end

%% TFR Low Frequency Plotting
% conditions = 'DifFB';
% proc_id   = 'eeg_full_ft';
% an_ids     = 'TFR_Fz_F2t1_z2t0_fl2t14';
% save_fig  = 1;
% fig_vis   = 'on';
% fig_ftype = 'png';
% 
% %% Compare p values across analyses
% proc_id = 'eeg_full_ft';
% an_ids   = 'ERP_Fz_F2t1_dm2t0_fl05t20';%'ERP_Pz_F2t1_dm2t0_fl05t20';
% do_id   = 'DifOut_lme_st0t5';
% do_lab  = {'Dif','Out','Dif*Out'};
% rl_id   = 'RL_DO_lme_st0t5';
% rl_lab  = {'pWin','sPE','uPE'};
% 
% % Load DO p values
% load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' do_id '_' an_ids '.mat']);
% do_pvals = nan([numel(do_lab) numel(lme)]);
% for grp_ix = 1:numel(do_lab)
%     for t_ix = 1:numel(lme)
%         do_pvals(grp_ix,t_ix) = lme{t_ix}.Coefficients.pValue(grp_ix+1);
%     end
% end
% 
% % Load RL p values
% load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' rl_id '_' an_ids '.mat']);
% rl_pvals = nan([numel(rl_lab) numel(lme)]);
% for grp_ix = 1:numel(rl_lab)
%     for t_ix = 1:numel(lme)
%         rl_pvals(grp_ix,t_ix) = lme{t_ix}.Coefficients.pValue(grp_ix+1);
%     end
% end
% 
% % Load time vector
% eval(['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' do_id '_vars.m']);
% eval(['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m']);
% load([SBJ_vars.dirs.proc,SBJs{1},'_',an_ids,'.mat']);
% cfgs = []; cfgs.latency = st.stat_lim;
% roi = ft_selectdata(cfgs, roi);
% time_vec = roi.time{1};
% 
% % Plot p value comparison
% figure;
% for grp_ix = 1:3
%     subplot(2,3,grp_ix); hold on;
%     plot(time_vec,do_pvals(grp_ix,:),'color','r');
%     plot(time_vec,rl_pvals(grp_ix,:),'color','k');
%     legend(do_lab{grp_ix},rl_lab{grp_ix});
% end
% for grp_ix = 1:3
%     subplot(2,3,3+grp_ix);
%     plot(time_vec,do_pvals(grp_ix,:)-rl_pvals(grp_ix,:),'color','b');
%     legend('DO - RL');
%     ylabel('DO better ... RL better');
% end