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
SBJs = {'EEG01','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};%,
% Bad SBJs:
%   'EEG02' (no prototype)
%   ,'EEG05' (no RLpRTlD model built!)

%% Prototype Selection
proc_id   = 'odd_full_ft';
cpa_id    = 'CPA';
plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 1: numel(SBJs)
    %ODD02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 'odd_full_ft', 1, 'on', 'ERP_stack_full_events_odd')
    %SBJ02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 1, 'on')
    %SBJ01_preproc(SBJs{x}, 'eeg_full_ft');
    SBJ06a_CPA_prototype_selection(SBJs{s}, proc_id, cpa_id, plt_id,...
        save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Test Candidate
eeg_proc_id = 'eeg_full_ft';
odd_proc_id = 'odd_full_ft';
conditions  = 'DifFB';
cpa_id      = 'CPA';
%an_ids      = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};%{'ERP_Fz_S15t28_dm2t0_fl05t20','ERP_Pz_S15t28_dm2t0_fl05t20'};
an_ids      = {'POW_Fz_F2t1_dm2t0_fl4t8','POW_Pz_F2t1_dm2t0_fl1t3'};
stat_id     = 'RLpRTlD_all_glm_st0t5';
save_fig    = 1;
fig_vis     = 'on';
fig_ftype   = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
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
%         % Run RL stats on single subject condidates
%         SBJ06d_CPA_candidate_stats_RL_SBJ(SBJs{s},eeg_proc_id,cpa_id,an_ids{an_ix},stat_id);
        
        % Plot RL Model Results
        plt_id    = 'ts_F2to1_evnts_sigLine';
        SBJ06e_CPA_candidate_ERP_plot_RL_fits(SBJs{s},eeg_proc_id,cpa_id,an_ids{an_ix},stat_id,...
            plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    close all
end

