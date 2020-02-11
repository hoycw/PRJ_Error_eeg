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
SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
           'EEG01','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12'};%'EEG02',
tfr_SBJs = {'EP06','EP07','EP08','EP10','EP16',...%,'EP11','EP14','EP15','EP17','EP18','EP19'
           'EEG01','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12'};%'EEG02',
% Bad SBJ:
%   EP01, EP02, EP05- recording errors
%   EP03- low quality
%   EP04- weird behavior?
%   EP09- 2 BDFs, unknown quality?
%   EP12, 13- don't exist
%   EP15- low quality?
%   EEG02- low quality
%   EEG09- multiple blocks, needs redo???
%   EEG11- recording failure

%% Single SBJ RL Model
proc_id   = 'eeg_full_ft';
stat_ids  = {'RLpRTulD_all_lme_st0t5'};%,'RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5','RLpRTlD_all_lme_st0t5'};
% RL models:
%   RL/pWinPEus (original) = pWin, sPE, uPE
%   RLRT = RL + tRT
%   RLpRT = RLRT + ptRT
%   RLpRTlD = RLpRT + lDist
%   RLpRTiD = RLpRT + iDist

for s = 1:numel(SBJs)
    for st_ix = 1:numel(stat_ids)
        SBJ04a_RL_model(SBJs{s},proc_id,stat_ids{st_ix});
    end
    close all;
end

%% ERP: Linear Mixed Effects Model (Over Time)
% % Main RL Model
% an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
% stat_ids  = {'RLpRTulD_all_lme_st0t5'};%'RLpRTulD_all_lme_st0t5','RLpRTlD_all_lme_st0t5','RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5',
% plt_id    = 'ts_F2to1_evnts_sigLine';
% Pre-Feedback RL Model
an_ids    = {'ERP_Fz_F4t1_dm4t3_fl05t20','ERP_Pz_F4t1_dm4t3_fl05t20'};
stat_ids  = {'RLpRTulD_all_lme_st3t5'};%'RLpRTulD_all_lme_st0t5','RLpRTlD_all_lme_st0t5','RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5',
plt_id    = 'ts_F4t1_evnts_sigLine';%'ts_F2to1_evnts_sigLine';

proc_id   = 'eeg_full_ft';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
%         SBJ04c_ERP_grp_stats_LME_RL(SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ04d_ERP_plot_stats_LME_RL_fits(SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
%     % Model Comparison Plots (Adjusted R-Squared)
%     SBJ04e_ERP_plot_RL_model_comparison(SBJs,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% ERP: Linear Mixed Effects Model (Mean Windows)
proc_id   = 'eeg_full_ft';
an_ids    = {'ERP_all_F2t1_dm2t0_fl05t20'};
stat_ids  = {'RLpRTlD_all_lme_mn2t3','RLpRTlD_all_lme_mn3t4'};
plt_ids   = {'topo_F18t25','topo_F3t45'};
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        SBJ04c_ERP_grp_stats_LME_RL(SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix});
%         SBJ04d_ERP_plot_stats_LME_RL_topo_reg(SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix},...
%             plt_ids{st_ix},save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Model Comparison Plots (Adjusted R-Squared)
%     error('need topo version!');
%     SBJ04e_ERP_plot_RL_model_comparison(SBJs,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% TFR: Linear Mixed Effects Model (Over Time)
% conditions = 'DifFB';
proc_id   = 'eeg_full_ft';
an_ids     = {'TFR_Fz_F2t1_db2t0_fl1t12b05','TFR_Pz_F2t1_db2t0_fl1t12b05'};
stat_ids  = {'RLpRTulD_all_lme_st0t5'};%'RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5',
%plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        SBJ05d_TFR_grp_stats_LME_RL(tfr_SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ05e_TFR_plot_stats_LME_RL_fits(tfr_SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix},save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Model Comparison Plots (Adjusted R-Squared)
%     error('this needs to be a martix version!');
%     SBJ05f_TFR_plot_RL_model_comparison(SBJs,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% POW: Linear Mixed Effects Model (Over Time)
proc_id   = 'eeg_full_ft';
an_ids     = {'POW_Fz_F2t1_db2t0_fl4t8','POW_Fz_F2t1_db2t0_fl8t12','POW_Pz_F2t1_db2t0_fl1t4'};
%topo_plt_ids = {'topo_F18t25','topo_F18t25','topo_F3t45'};
stat_ids  = {'RLpRTulD_all_lme_st0t5'};%'RL_all_lme_st0t5','RLRT_all_lme_st0t5','RLpRT_all_lme_st0t5',
plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_ids)
        SBJ05d_TFR_grp_stats_LME_RL(SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ05e_POW_plot_stats_LME_RL_fits(SBJs,proc_id,an_ids{an_ix},stat_ids{st_ix},plt_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Model Comparison Plots (Adjusted R-Squared)
%     SBJ05f_POW_plot_RL_model_comparison(SBJs,proc_id,an_ids{an_ix},stat_ids,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

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