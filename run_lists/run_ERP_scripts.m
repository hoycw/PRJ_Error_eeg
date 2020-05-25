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
SBJ_id = 'goodall';
SBJs = fn_load_SBJ_list(SBJ_id);

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

%% ERPs: Fz and Pz
%   RL Model Analysis:
an_ids     = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
conditions = 'EHSu';
proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';

for an_ix = 1%:numel(an_ids)
    plt_id     = 'stack_F2t1_evnt_c5';
    for s = 1:numel(SBJs)
%         SBJ03a_ERP_save(SBJs{s},proc_id,an_ids{an_ix});
%         SBJ03b_ERP_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%               'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         SBJ03b_ERP_plot_butterfly(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         SBJ03b_ERP_plot_stack(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    SBJ03c_ERP_plot_grp(SBJ_id,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
%     plt_id = 'ts_F2to1_but_evnts_sigPatch';
%     SBJ03c_ERP_plot_grp_butterfly(SBJs,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     close all;
end

%% Plot only surprise ERPs
an_ids     = {'ERP_Fz_F2t1_dm2t0_fl05t20'};
conditions = 'EHSu';
proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';
plt_id     = 'ts_F1t5_evnts_sigLine';

for an_ix = 1%:numel(an_ids)
    SBJ03c_ERP_plot_grp(SBJ_id,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% ERPs: Full Cap Topography
an_id      = 'ERP_all_S2t1_dm2t0_fl05t20';

proc_id    = 'eeg_full_ft';
conditions = 'DifFB';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 1:numel(SBJs)
    SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
%     % FRN by condition
%     plt_id    = 'topo_F18t25';
%     SBJ03b_ERP_plot_topo_cond(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     
%     % P3 by condition
%     plt_id    = 'topo_F3t45';
%     SBJ03b_ERP_plot_topo_cond(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     close all;
end

% % FRN Group plot
% plt_id    = 'topo_F18t25';
% SBJ03c_ERP_plot_grp_topo_cond(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% 
% % P3 Group Plot
% plt_id    = 'topo_F3t45';
% SBJ03c_ERP_plot_grp_topo_cond(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%% ========================================================================
%   OLD UNUSED ANALYSES (not going in the paper)
%  ========================================================================

% %% View difference wave ERPs
% % proc_id    = 'eeg_full_ft';
% % an_id      = 'ERP_Cz_F2t1_dm2t0_fl05t20';
% % conditions = 'DifOutS';
% % plt_id     = 'ts_F2to1_evnts_sigLine';
% % save_fig   = 1;
% % fig_vis    = 'on';
% % fig_ftype  = 'png';
% % for s = 1:numel(SBJs)
% %     SBJ03b_ERP_plot_diffwave(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
% %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% % end
% % SBJ03c_ERP_plot_grp_diffwave(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
% %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% % 
% % plt_id = 'ts_F2to1_but_evnts_sigPatch';
% % SBJ03c_ERP_plot_grp_diffwave_butterfly(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
% %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);

