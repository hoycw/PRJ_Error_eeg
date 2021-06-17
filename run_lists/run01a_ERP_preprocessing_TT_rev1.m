%% ERP analysis and plotting  of Target Time task for Sequential PE Revision 1 Submission
% Developed over time, but editted 4/21/21 by Colin W Hoy
% 	Sup. Fig. 10A: SBJ03c_ERP_plot_grp_pkLine

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
SBJ_id = 'goodall';%'good1';%
SBJs = fn_load_SBJ_list(SBJ_id);

%% ERPs: Fz and Pz over time
%   RL Model Analysis:
an_ids     = {'ERP_Fz_F2t1_dm2t0_fl05t20'};%,'ERP_Pz_F2t1_dm2t0_fl05t20'};
conditions = 'DifWL';%'DifFB';%'EHSu';
proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';%'svg';%

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        % Re-align data to event, select channels and epoch, filter, save
        %   Options to downsample and run LaPlacian transform
%         SBJ03a_ERP_save(SBJs{s},proc_id,an_ids{an_ix});

        % Plot SBJ ERPs per condition
%         plt_id = 'ts_F2to1_evnts_sigLine';
%         SBJ03b_ERP_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%               'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot SBJ single trials, one subplot per condition
%         SBJ03b_ERP_plot_butterfly(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);

        % Plot SBJ single trial stack with ERP underneath
%         plt_id     = 'stack_F2t1_evnt_c5';
%         SBJ03b_ERP_plot_stack(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Plot Group ERPs
    %   can reproduce plots in Fig. 2A (Fz) and 2B (Pz), but see SBJ04d_ERP_plot_stats_LME_RL_fits.m
    SBJ03c_ERP_plot_grp(SBJ_id,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Plot all SBJ ERPs overlapping (butterfly)
%     plt_id = 'ts_F2to1_but_evnts_sigPatch';
%     SBJ03c_ERP_plot_grp_butterfly(SBJs,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     close all;
end

%% Plot ERPs with FRN peaks marked
%*** plots Sup. Fig. 10A (easy/hard FB separately, neutral for Easy+Hard)
an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';
cond_list  = {'HdOutS','EzOutS','EHSu'};
proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';
plt_id     = 'ts_F0t5_evnts_sigLine';

% FRN Peak Detection Settings
pk_lim  = [0.18 0.3];
pk_sign = -1;

for cond_ix = 1:numel(cond_list)
    SBJ03c_ERP_plot_grp_pkLine(SBJ_id,cond_list{cond_ix},proc_id,an_id,pk_lim,pk_sign,...
        plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Save Group ERPs to get data-driven analysis windows from ERP peak times
% NOTE: Group FRN peak time across all conditions/subjects is used to
%   center traditional mean window analyses using:
%       -SBJ04c_ERP_grp_stats_LME_mean_window
%       -SBJ04d_ERP_plot_stats_LME_mean_betas
conditions = 'All';
an_id      = 'ERP_Fz_F2t1_dm2t0_fl05t20';%'ERP_Pz_F2t1_dm2t0_fl05t20';
proc_id    = 'eeg_full_ft';

SBJ03c_ERP_save_grp_ERP_cond(SBJ_id,conditions,proc_id,an_id);

% NOTE: SBJ03c_ERP_save_grp_topo_cond does something similar but seems to
% be unused... maybe it was a precursor to SBJ03c_ERP_plot_grp_topo_ts_cond?

%% ERPs: Full Cap Topography
% NOTE: Sup. Fig. 6 (ERP topo dynamics) plotted in run03_RL_model_results
% because it locks the averaging windows to peaks of model coefficients
an_id      = 'ERP_all_F2t1_dm2t0_fl05t20';

proc_id    = 'eeg_full_ft';
conditions = 'DifFB';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 1:numel(SBJs)
    % Compute ERP across all channels
    SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
    
%     % Plot FRN window topo by condition
%     plt_id    = 'topo_F18t25';
%     SBJ03b_ERP_plot_topo_cond(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     
%     % Plot P3 window topo by condition
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


