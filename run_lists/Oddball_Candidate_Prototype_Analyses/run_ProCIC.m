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
<<<<<<< HEAD
SBJ_id = 'goodall';
SBJs = fn_load_SBJ_list(SBJ_id);
=======
SBJ_id = 'good2';%'goodEEG';
sbj_file = fopen([root_dir 'PRJ_Error_EEG/scripts/SBJ_lists/' SBJ_id '.sbj']);
tmp = textscan(sbj_file,'%s');
fclose(sbj_file);
SBJs = tmp{1}; clear tmp;
>>>>>>> fa78d68822e09830587ce6d033850dc2563ba55c

%% Generate condition-specific topos
conditions = 'Odd';
proc_id    = 'odd_full_ft';
an_id      = 'ERP_all_S2t1_dm2t0_fl05t20';
SBJ03c_ERP_save_grp_topo_cond(SBJ_id,conditions,proc_id,an_id);

%% Prototype Selection
proc_id   = 'odd_full_ft';
cpa_id    = 'CPA_main';
an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 1:numel(SBJs)
    %ODD02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 'odd_full_ft', 1, 'on', 'ERP_stack_full_events_odd')
    %SBJ02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 1, 'on')
    %SBJ01_preproc(SBJs{x}, 'eeg_full_ft');
    %     SBJ06a_CPA_prototype_selection(SBJs{s}, proc_id, 'CPA', plt_id,...
    %         save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Topo correlation
    SBJ06a_CPA_prototype_selection(SBJs{s}, proc_id, cpa_id, save_fig);
    SBJ06b_CPA_prototype_plot(SBJs{s}, proc_id, cpa_id, plt_id,...
        save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    close all;
end

%% Reconstruct Prototype Oddball ERPs
%-------------- Time series --------------
proc_id    = 'odd_full_ft';
conditions = 'Odd';
cpa_id     = 'CPA_main';
an_ids     = {'ERP_Fz_S2t1_dm2t0_fl05t20','ERP_Pz_S2t1_dm2t0_fl05t20'};
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        SBJ06b_CPA_prototype_ERP_save(SBJs{s}, proc_id, an_ids{an_ix}, cpa_id);
    end
    plt_id    = 'ts_S2t1_evnts_sigLine';
    SBJ06c_CPA_prototype_ERP_plot_GRP(SBJ_id,conditions,proc_id,cpa_id,an_ids{an_ix},plt_id,...
        save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%-------------- Full Cap --------------
topo_an_id  = 'ERP_all_S2t1_dm2t0_fl05t20';
for s = 1:numel(SBJs)
    SBJ06b_CPA_prototype_ERP_save(SBJs{s}, proc_id, topo_an_id, cpa_id);
end
plt_id = 'topo_F3t45';
SBJ06c_CPA_prototype_ERP_plot_GRP_topo_cond(SBJ_id,conditions,proc_id,cpa_id,topo_an_id,plt_id,...
    save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%% Print # prototypes per SBJ
for s = 1:numel(SBJs)
%     load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_' cpa_id '_' proc_id '_prototype.mat']);
%     fprintf('%s: %d\n',SBJs{s},numel(final_ics));
end

%% Test Candidate
eeg_proc_id = 'eeg_full_ft';
odd_proc_id = 'odd_full_ft';
conditions  = 'DifFB';
cpa_id      = 'CPA_main';
an_ids      = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
stat_id     = 'RL_all_lme_st05t5';
save_fig    = 1;
fig_vis     = 'on';
fig_ftype   = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        % Compute ERP for Candidate
%         SBJ07a_CPA_candidate_ERP_save(SBJs{s},eeg_proc_id,odd_proc_id,an_ids{an_ix},cpa_id);
        
        % Plot full TT epoch with trial stack
        plt_id    = 'stack_F2t1_evnt_c5';%'stack_S15to28_evnt_c5';
%         SBJ07b_CPA_candidate_ERP_plot_stack(SBJs{s},conditions,eeg_proc_id,cpa_id,an_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot feedback TT epoch with only ERPs
%         plt_id    = 'ts_F2to1_evnts_sigLine';
%         SBJ07b_CPA_candidate_ERP_plot(SBJs{s},conditions,eeg_proc_id,cpa_id,an_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         close all;
    end
    
    % Group LME Stats
%     SBJ07c_CPA_candidate_ERP_GRP_stats_LME_RL(SBJ_id,eeg_proc_id,cpa_id,an_ids{an_ix},stat_id);
    
    plt_id    = 'ts_F2to1_evnts_sigLine';
    SBJ07d_CPA_candidate_ERP_plot_LME_RL_fits(SBJ_id,eeg_proc_id,cpa_id,an_ids{an_ix},stat_id,...
        plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    %close all
end

%% Topo ERPs:
eeg_proc_id = 'eeg_full_ft';
odd_proc_id = 'odd_full_ft';
conditions  = 'DifFB';
cpa_id      = 'CPA_main';
an_id       = 'ERP_all_F2t1_dm2t0_fl05t20';
stat_ids    = {'RL_all_lme_mn1sPE','RL_all_lme_mn1uPE'};
plt_ids     = {'topo_F18t25','topo_F3t45'};
save_fig    = 1;
fig_vis     = 'on';
fig_ftype   = 'png';

for st_ix = 1:numel(stat_ids)
    for s = 1:numel(SBJs)
        % Compute ERP for Candidate
        SBJ07a_CPA_candidate_ERP_save(SBJs{s},eeg_proc_id,odd_proc_id,an_id,cpa_id);
    end
    
    % Group LME Stats
    SBJ07c_CPA_candidate_ERP_GRP_stats_LME_RL(SBJ_id,eeg_proc_id,cpa_id,an_id,stat_ids{st_ix});
    
    SBJ07d_CPA_candidate_ERP_plot_LME_RL_topo_reg(SBJ_id,cpa_id,an_id,stat_ids{st_ix},...
        plt_ids{st_ix},save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    %close all
end

%% Run ERP plots for SBJ by SBJ comparison
% proc_id   = 'eeg_full_ft';
% conditions = 'DifFB';
% an_ids      = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
% plt_id    = 'ts_F2to1_evnts_sigLine';
% save_fig    = 1;
% fig_vis     = 'on';
% fig_ftype   = 'png';
% 
% for s = 1:numel(SBJs)
% %     SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
%     SBJ03b_ERP_plot(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% end

%% Power TFR: Linear Mixed Effects Model (Over Time-Frequency Power)
eeg_proc_id = 'eeg_full_ft';
odd_proc_id = 'odd_full_ft';
cpa_id      = 'CPA_main';
an_ids      = {'TFR_Fz_F2t1_db2t0_fl1t12','TFR_Pz_F2t1_db2t0_fl1t12'};
stat_ids    = {'RL_all_lme_st0t5'};
save_fig    = 1;
fig_vis     = 'on';
fig_ftype   = 'png';

for an_ix = 1:numel(an_ids)
    % Compute TFRs
    for s = 1:numel(SBJs)
%         SBJ08a_CPA_candidate_TFR_save(SBJs{s}, eeg_proc_id, odd_proc_id, cpa_id, an_ids{an_ix});
    end
    
    for st_ix = 1:numel(stat_ids)
        SBJ08d_CPA_candidate_TFR_grp_stats_LME_RL(SBJ_id,eeg_proc_id,cpa_id,an_ids{an_ix},stat_ids{st_ix});
        SBJ08e_CPA_candidate_TFR_plot_stats_LME_RL_fits(SBJ_id,eeg_proc_id,cpa_id,an_ids{an_ix},stat_ids{st_ix},save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
end

%% Phase TFR: Linear Mixed Effects Model (Over Time-Frequency Phase)
eeg_proc_id = 'eeg_full_ft';
odd_proc_id = 'odd_full_ft';
cpa_id      = 'CPA_main';
an_ids      = {'PHS_Fz_F2t1_fl1t12'};
stat_ids    = {'RL_all_lme_st0t5'};
save_fig    = 1;
fig_vis     = 'on';
fig_ftype   = 'png';

for an_ix = 1:numel(an_ids)
    % Compute TFRs
    for s = 1:numel(SBJs)
        SBJ08a_CPA_candidate_TFR_save(SBJs{s}, eeg_proc_id, odd_proc_id, cpa_id, an_ids{an_ix});
    end
    
    for st_ix = 1:numel(stat_ids)
%         SBJ08d_CPA_candidate_TFR_grp_stats_LME_RL(SBJ_id,eeg_proc_id,cpa_id,an_ids{an_ix},stat_ids{st_ix});
%         SBJ08e_CPA_candidate_TFR_plot_stats_LME_RL_fits(SBJ_id,eeg_proc_id,cpa_id,an_ids{an_ix},stat_ids{st_ix},save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
end


