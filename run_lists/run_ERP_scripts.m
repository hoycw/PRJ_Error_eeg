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
SBJs = {'EP07','EP08','EP10','EP11','EP14','EP16','EP17','EP19',...
           'EEG01','EEG03','EEG04','EEG05','EEG06','EEG08','EEG10'};
% Not Ready SBJ:
%   EP06: only 62 channels?
%   EP09: 2 BDFs, unknown quality?
%   EP15: low quality?
%   EP18: low trial count (328)
%   EEG07: low trial count (373)
%   EEG09: multiple blocks, needs redo???
%   EEG12: low trial count (271)
%   EEG13-27: ready, not used yet
% Bad SBJ:
%   EP01, EP02, EP05: recording errors
%   EP03: low quality
%   EP04: weird behavior?
%   EP12-13: don't exist
%   EEG02: low quality
%   EEG11: recording failure

%% Run preprocessing
proc_id_ica = 'eeg_full_ft';
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
%   RL Model Analysis:
an_ids     = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
conditions = 'DifFB';
% %   Pre-Feedback RL Model Analysis:
% an_ids     = {'ERP_Fz_F4t1_dm4t3_fl05t20','ERP_Pz_F4t1_dm4t3_fl05t20'};
% conditions = 'DifFB';
%   QA Plotting:
% an_id      = 'ERP_Z4_F2t1_dm2t0_fl05t20';
% conditions = 'Dif';%FB

proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    plt_id     = 'ts_F2to1_evnts_sigLine';
    for s = 1:numel(SBJs)
        SBJ03a_ERP_save(SBJs{s},proc_id,an_ids{an_ix});
%         SBJ03b_ERP_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
%     SBJ03c_ERP_plot_grp(SBJs,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
%     plt_id = 'ts_F2to1_but_evnts_sigPatch';
%     SBJ03c_ERP_plot_grp_butterfly(SBJs,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    close all;
end

%% Plot ERP Topos
proc_id    = 'eeg_full_ft';
conditions = 'DifFB';
an_id      = 'ERP_all_F2t1_dm2t0_fl05t20';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 1:numel(SBJs)
    % FRN by condition
    plt_id    = 'topo_F18t25';
    SBJ03b_ERP_plot_topo_cond(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % P3 by condition
    plt_id    = 'topo_F3t45';
    SBJ03b_ERP_plot_topo_cond(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    close all;
end

% FRN Group plot
plt_id    = 'topo_F18t25';
SBJ03c_ERP_plot_grp_topo_cond(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

% P3 Group Plot
plt_id    = 'topo_F3t45';
SBJ03c_ERP_plot_grp_topo_cond(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%% View difference wave ERPs
% proc_id    = 'eeg_full_ft';
% an_id      = 'ERP_Cz_F2t1_dm2t0_fl05t20';
% conditions = 'DifOutS';
% plt_id     = 'ts_F2to1_evnts_sigLine';
% save_fig   = 1;
% fig_vis    = 'on';
% fig_ftype  = 'png';
% for s = 1:numel(SBJs)
%     SBJ03b_ERP_plot_diffwave(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% end
% SBJ03c_ERP_plot_grp_diffwave(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% 
% plt_id = 'ts_F2to1_but_evnts_sigPatch';
% SBJ03c_ERP_plot_grp_diffwave_butterfly(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%% ERP Stats: Window Mean
proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';
% for s = 1:numel(SBJs)
%     SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
% end

an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_id    = 'DifOut_anv_mn2t3_jk';
SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_id    = 'DifFB_anv_mn2t3_jk';
SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

an_id      = 'ERP_Pz_F2t1_dm2t0_fl05t20';
stat_id    = 'DifOut_anv_mn3t4_jk';
SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

an_id      = 'ERP_Pz_F2t1_dm2t0_fl05t20';
stat_id    = 'DifFB_anv_mn3t4_jk';
SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%% ERP Stats: Peak-to-Peak
proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

an_id   = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_id = 'DifOut_anv_p2pFRN';
SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

an_id   = 'ERP_Fz_F2t1_dm2t0_fl05t20';
stat_id = 'DifFB_anv_p2pFRN';
SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);

