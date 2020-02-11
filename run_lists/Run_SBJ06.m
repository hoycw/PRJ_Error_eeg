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
SBJs = {'EEG01','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG12'};%'EEG03','EEG10',
% Bad SBJs:
%   EPs: no oddball
%   'EEG02': bad quality

%% Prototype Selection
proc_id   = 'odd_full_ft';
% cpa_id    = 'CPA';%'CPA_odd';
plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 2:numel(SBJs)
    %ODD02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 'odd_full_ft', 1, 'on', 'ERP_stack_full_events_odd')
    %SBJ02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 1, 'on')
    %SBJ01_preproc(SBJs{x}, 'eeg_full_ft');
%     SBJ06a_CPA_prototype_selection(SBJs{s}, proc_id, 'CPA', plt_id,...
%         save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    try
        SBJ06a_CPA_prototype_selection(SBJs{s}, proc_id, 'CPA_odd2', plt_id,...
            save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);

    catch
        fprintf(2,'No ICs found for %s\n',SBJs{s});
    end
end

for s = 1:numel(SBJs)
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_CPA_' proc_id '_prototype.mat']);
    orig_ics = final_ics;
    try
        load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_CPA_odd2_' proc_id '_prototype.mat']);
    catch
        final_ics = [];
    end
    fprintf('%s: n_tar = %d; n_odd2 = %d; overlap = %d\n',SBJs{s},numel(orig_ics),numel(final_ics),numel(intersect(orig_ics,final_ics)));
end
% COMPARISON OUTPUT:
% EEG01: n_origs = 3; n_odd = 3; overlap = 3
% EEG03: n_origs = 1; n_odd = 1; overlap = 0
% EEG04: n_origs = 1; n_odd = 1; overlap = 1
% EEG05: n_origs = 2; n_odd = 1; overlap = 0
% EEG06: n_origs = 2; n_odd = 4; overlap = 2
% EEG07: n_origs = 1; n_odd = 2; overlap = 1
% EEG08: n_origs = 1; n_odd = 2; overlap = 1
% EEG09: n_origs = 4; n_odd = 6; overlap = 3
% EEG10: n_origs = 1; n_odd = 1; overlap = 0
% EEG12: n_origs = 2; n_odd = 2; overlap = 2
% 
%========
% EEG01: n_origs = 1; n_odd = 0; overlap = 0
% EEG03: n_origs = 1; n_odd = 1 (but #42); overlap = 1
% EEG04: n_origs = 1 (but #25); n_odd = 1 (but same #25); overlap = 1
% EEG05: n_origs = 1; n_odd = 0; overlap = 0
% EEG06: n_origs = 1; n_odd = 1; overlap = 1
% EEG07: n_origs = 1; n_odd = 0; overlap = 0
% EEG08: n_origs = 1; n_odd = 2; overlap = 1
% EEG09: n_origs = 1; n_odd = 1; overlap = 1
% EEG10: n_origs = 1; n_odd = 1 (but #40); overlap = 1
% EEG12: n_origs = 1; n_odd = 1; overlap = 1

%% Test Candidate
eeg_proc_id = 'eeg_full_ft';
odd_proc_id = 'odd_full_ft';
conditions  = 'DifFB';
cpa_id      = 'CPA';
an_ids      = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};%{'ERP_Fz_S15t28_dm2t0_fl05t20','ERP_Pz_S15t28_dm2t0_fl05t20'};
% an_ids      = {'POW_Fz_F2t1_dm2t0_fl4t8','POW_Pz_F2t1_dm2t0_fl1t3'};
stat_id     = 'RLpRTulD_all_lme_st0t5';
save_fig    = 1;
fig_vis     = 'on';
fig_ftype   = 'png';

for an_ix = 1:numel(an_ids)
%     for s = 1:numel(SBJs)
%         % Compute ERP for Candidate
%         SBJ06b_CPA_candidate_ERP_save(SBJs{s},eeg_proc_id,odd_proc_id,an_ids{an_ix},cpa_id);
%         
%         % Plot full TT epoch with trial stack
%         plt_id    = 'stack_F2t1_evnt_c5';%'stack_S15to28_evnt_c5';
%         SBJ06c_CPA_candidate_ERP_plot_stack(SBJs{s},conditions,eeg_proc_id,cpa_id,an_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         
%         % Plot feedback TT epoch with only ERPs
%         plt_id    = 'ts_F2to1_evnts_sigLine';
%         SBJ06c_CPA_candidate_ERP_plot(SBJs{s},conditions,eeg_proc_id,cpa_id,an_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         
% %         % Run RL stats on single subject condidates
% %         SBJ06d_CPA_candidate_stats_RL_SBJ(SBJs{s},eeg_proc_id,cpa_id,an_ids{an_ix},stat_id);
% %         
% %         % Plot RL Model Results
% %         plt_id    = 'ts_F2to1_evnts_sigLine';
% %         SBJ06e_CPA_candidate_ERP_plot_RL_fits(SBJs{s},eeg_proc_id,cpa_id,an_ids{an_ix},stat_id,...
% %             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         close all;
%     end
    
    % Group LME Stats
%     SBJ06d_CPA_candidate_ERP_GRP_stats_LME_RL(SBJs,eeg_proc_id,cpa_id,an_ids{an_ix},stat_id);
    
    plt_id    = 'ts_F2to1_evnts_sigLine';
    SBJ06e_CPA_candidate_ERP_plot_LME_RL_fits(SBJs,eeg_proc_id,cpa_id,an_ids{an_ix},stat_id,...
        plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
%     % Group plotting of single SBJ stats
%     plt_id    = 'ts_F0t5_but_evnts_sigPatch';
%     SBJ06f_CPA_candidate_ERP_plot_GRP_betas(SBJs,eeg_proc_id,cpa_id,an_ids{an_ix},stat_id,...
%         plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    %close all
end

%% Run ERP plots for SBJ by SBJ comparison
proc_id   = 'eeg_full_ft';
conditions = 'DifFB';
an_ids      = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};%{'ERP_Fz_S15t28_dm2t0_fl05t20','ERP_Pz_S15t28_dm2t0_fl05t20'};
%an_id     = 'POW_FCz_F2t1_dm2t0_fl4t8';%'ERP_FCz_F2t1_dm2t0_fl05t20';%'POW_all_F2t1_dm2t0_fl4t8';%'ERP_all_F2t1_dm2t0_fl05t20';
plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig    = 1;
fig_vis     = 'on';
fig_ftype   = 'png';

for s = 1:numel(SBJs)
%     SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
    SBJ03b_ERP_plot(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

