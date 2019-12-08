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
           'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};

%% Plot Accuracy
cond_lab = {'easy','hard'};
acc_cond = zeros([numel(SBJs) 2]);
for s = 1:numel(SBJs)
    [bhv] = fn_load_behav_csv([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' SBJs{s} '_behav.csv']);
    for cond_ix = 1:2
        acc_cond(s,cond_ix) = mean(bhv.hit(strcmp(bhv.cond,cond_lab{cond_ix})));
    end
end
fprintf('Accuracy (n=%i):\n',numel(SBJs));
for cond_ix = 1:2
    fprintf('\t%s = %f +/- %f\n',cond_lab{cond_ix},mean(acc_cond(:,cond_ix)),std(acc_cond(:,cond_ix)));
end

%% Run preprocessing
proc_id_ica = proc_id;
gen_figs    = 0;
fig_vis     = 'off';
reject_visual = 0;
plot_final_check = 0;

% SBJ_times = zeros(size(SBJs));
% tic;
% for s = 1:numel(SBJs)
%     SBJ02a_artifact_rejection(SBJs{s}, proc_id, gen_figs, fig_vis)
%     SBJ02b_ica_rejection(SBJs{s}, proc_id, proc_id_ica, reject_visual);
%     SBJ02c_trial_rejection(SBJs{s}, proc_id, plot_final_check)
%     SBJ_times(s) = toc;
%     if s==1; elapsed = SBJ_times(s); else; elapsed = SBJ_times(s)-SBJ_times(s-1); end
%     fprintf('%s preprocessed at %.1f s (SBJ time = %.1f)\n',SBJs{s},SBJ_times(s),elapsed);
% end

%% View basic ERPs
proc_id    = 'eeg_full_ft';
an_ids     = {'ERP_Fz_F2t1_dm2t0_fl05t20'};%,'ERP_Pz_F2t1_dm2t0_fl05t20'
stat_conds = {'EzOutS','HdOutS'};%'DifOut','DifFB'};%,'DifOutS'};
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';

for an_ix = 1:numel(an_ids)
    for st_ix = 1:numel(stat_conds)
        plt_id     = 'ts_F2to1_evnts_sigLine';
        %     for s = 1:numel(SBJs)
        %         SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
        %         SBJ03b_ERP_plot(SBJs{s},stat_conds{st_ix},proc_id,an_id,plt_id,save_fig,...
        %             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        %     end
        SBJ03c_ERP_plot_grp(SBJs,stat_conds{st_ix},proc_id,an_ids{an_ix},plt_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        %     plt_id = 'ts_F2to1_but_evnts_sigPatch';
        %     SBJ03c_ERP_plot_grp_butterfly(SBJs,stat_conds{st_ix},proc_id,an_id,plt_id,save_fig,...
        %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        %
        %     close all;
    end
end

%% View difference wave ERPs
% proc_id    = 'eeg_full_ft';
% an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
% stat_conds = {'DifOutdO','DifOutWL','DifOutUE'};
% plt_id     = 'ts_F2to1_evnts_sigLine';
% save_fig   = 1;
% fig_vis    = 'on';
% fig_ftype  = 'svg';
% % for s = 1:numel(SBJs)
% %     SBJ03b_ERP_plot_diffwave(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
% %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% % end
% for st_ix = 1:numel(stat_conds)
%     SBJ03c_ERP_plot_grp_diffwave(SBJs,stat_conds{st_ix},proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% end
% % plt_id = 'ts_F2to1_but_evnts_sigPatch';
% % SBJ03c_ERP_plot_grp_diffwave_butterfly(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
% %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%% ERP Stats: FRN Window Mean
proc_id    = 'eeg_full_ft';
an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_conds = {'DifOut_anv_mn2t3','DifOut_anv_mn2t3_jk','DifFB_anv_mn2t3_jk'};
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';

for st_ix = 1:numel(stat_conds)
    SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_conds{st_ix},save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% ERP Stats: P300 Window Mean
proc_id    = 'eeg_full_ft';
an_id      = 'ERP_Pz_F2t1_dm2t0_fl05t20';
stat_conds = {'DifOut_anv_mn3t4','DifOut_anv_mn3t4_jk',...
              'DifFB_anv_mn3t4','DifFB_anv_mn3t4_jk'};
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';

for st_ix = 1:numel(stat_conds)
    SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_conds{st_ix},save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% ERP Stats: FRN Peak-to-Peak
proc_id    = 'eeg_full_ft';
an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_conds = {'DifOut_anv_p2pFRN','DifFB_anv_p2pFRN'};
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';

for st_ix = 1:numel(stat_conds)
    SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_conds{st_ix},save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% FRN Peak latencies
proc_id    = 'eeg_full_ft';
an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_id    = 'DifFB_anv_p2pFRN';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';

SBJ04c_ERP_grp_stats_peak_times(SBJs,proc_id,an_id,stat_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
%% Plot TFRs
% conditions = 'DifOut';
% proc_id  = 'eeg_full_ft';
% an_id    = 'TFR_Fz_F2t1_rc2t0_fl2t40';
% save_fig = 1;
% 
% for s = 1:numel(SBJs)
%      SBJ05a_TFR_save(SBJs{s}, proc_id, an_id)
% %      SBJ05b_TFR_plot(SBJs{s}, 'DifOut', proc_id, an_id, save_fig);
% end
% 
% SBJ05c_TFR_plot_grp(SBJs,conditions,proc_id,an_id,save_fig);

%% Linear Mixed Effects Model
proc_id   = 'eeg_full_ft';
an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
stat_id   = 'DifOut_lme_st0t5';
plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

% for s = 1:numel(SBJs)
%     SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
% end

for an_ix = 1:numel(an_ids)
    % SBJ04c_ERP_grp_stats_LME(SBJs,proc_id,an_ids{an_ix},stat_id);
    SBJ04d_ERP_plot_stats_LME(SBJs,proc_id,an_ids{an_ix},stat_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    % SBJ04d_ERP_plot_stats_LME(SBJs,proc_id,an_ids{an_ix},stat_id,plt_id,save_fig,...
    %         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_median',1);
end

%% ODDBALL 
% odd_plt_id = 'ts_S2to13_evnts_sigPatch';
% an_id = 'ERP_Cz_F_trl15t28_flt05t20_stat06';
% odd_an_id = 'ERP_Cz_S2t13_flt5t20_st6';
% odd_proc_id = 'odd_full_ft';
% proc_id = 'eeg_full_ft';
% for s = 1:numel(SBJs)
%     PLOTSTAT_HitErr_Individ(SBJs{s}, proc_id, plt_id, an_id, 'off', 1, 'png');
%     PLOT_4Conds_Individ(SBJs{s}, proc_id, plt_id, an_id, 'off', 1, 'png');
%     if startsWith(SBJs{s},'EEG')
%         PLOTSTAT_Oddball_Individ(SBJs{s}, odd_proc_id, odd_plt_id, odd_an_id, 'off', 1, 'png');
%     end
%     %CALC_FRN_Individ(SBJs{s}, proc_id, plt_id, an_id);
% end
% 
% PLOTSTAT_HitErr_Group(SBJs, plt_id, an_id, 'off', 1, 'png');
% PLOT_4Conds_Group(SBJs,  plt_id, an_id, 'off', 1, 'png');
% PLOTSTAT_Oddball_Group(SBJs, odd_plt_id, odd_an_id, 'off', 1, 'png');
% %PLOTSTATS_FRN_Group(SBJs, proc_id, plt_id, an_id, 'off', 1, '.png');
% PLOT_DiffWave_Group(SBJs, plt_id, an_id, 'off', 1, '.png', 'DifWL');
% PLOT_DiffWave_Group(SBJs, plt_id, an_id, 'off', 1, '.png', 'DifEH');

%% Run Oddball ERPs
% odd_SBJs = {'EEG01','EEG02','EEG03','EEG04','EEG06','EEG08'};
% 
% proc_id   = 'odd_full_ft';
% an_id     = 'ERP_Cz_S2t13_flt5t20_st6';
% plt_id    = 'ts_S2to13_evnts_sigPatch';
% fig_vis   = 'on';
% save_fig  = 1;
% fig_ftype = 'png';
% 
% for s = 4:numel(SBJs)
% %     oddball_stats_plot(SBJs{s}, proc_id, plt_id, an_id, fig_vis, save_fig, fig_ftype)
% end