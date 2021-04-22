%% Time-frequency power and phase analyses for Target Time task in Sequential PE Revision 1 Submission
% Developed over time, but editted 4/16/21 by Colin W Hoy
% Oddball TFRs were not included in the revision 1 manuscript

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
SBJ_id = 'goodOB';%EEG1';%'goodall';%
SBJs = fn_load_SBJ_list(SBJ_id);

%% Compute TFRs for Power
ob_proc_id = 'odd_full_ft';
tt_proc_id = 'eeg_full_ft';
an_ids     = {'TFR_Fz_S2t1_db2t0_fl1t12', 'TFR_Pz_S2t1_db2t0_fl1t12'};
erp_ids    = {'ERP_Fz_S2t1_dm2t0_fl05t20','ERP_Pz_S2t1_dm2t0_fl05t20'};
conditions = 'OB';
plt_id     = 'ts_S2t8_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        % Reconstruct and clean raw data, filter, cut trials to event,
        %   baseline correct, select channels, save
%         SBJ07a_OB_TFR_save(SBJs{s}, ob_proc_id, tt_proc_id, an_ids{an_ix});
        
        % Plot TFRs of power data per condition for single SBJ
%         SBJ07b_OB_TFR_plot(SBJs{s}, conditions, ob_proc_id, an_ids{an_ix}, plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot TFRs of power data with ERP overlay per condition for single SBJ
%         SBJ07b_OB_TFR_ERP_plot(SBJs{s}, conditions, ob_proc_id, an_ids{an_ix},erp_ids{an_ix}, plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Plot TFRs of power data per condition for group
%     SBJ07c_OB_TFR_plot_grp(SBJ_id, conditions, ob_proc_id, an_ids{an_ix}, plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Plot TFRs of power data with ERP overlay per condition for group
%     SBJ07c_OB_TFR_ERP_plot_grp(SBJ_id, conditions, ob_proc_id, an_ids{an_ix},erp_ids{an_ix}, plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Compute Phase for ITPC
ob_proc_id = 'odd_full_ft';
tt_proc_id = 'eeg_full_ft';
an_ids     = {'PHS_Fz_S2t1_fl1t12'};%,'PHS_Pz_F2t1_fl1t12'};
erp_ids    = {'ERP_Fz_S2t1_dm2t0_fl05t20'};%,'ERP_Pz_F2t1_dm2t0_fl05t20'};
conditions = 'OB';
plt_id     = 'ts_S2t1_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        % Reconstruct and clean raw data, filter, cut trials to event,
        %   baseline correct, select channels, save
%         SBJ07a_OB_TFR_save(SBJs{s}, ob_proc_id, tt_proc_id, an_ids{an_ix});
        
        % Compute and plot ITPC with ERP overlay per condition for single SBJ
%         SBJ05b_ITC_ERP_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Exploratory rose plot: Compute and plot ITPC with ERP overlay per condition for single SBJ
        %   Also extracts/plots mean phase angle in T-F window across conditions
%         SBJ05b_ITC_ERP_rose_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},phs_freq_lim,phs_time_lim,erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Compute and plot ITPC with ERP overlay per condition for group
%     SBJ07c_OB_ITC_ERP_plot_grp(SBJ_id,conditions,ob_proc_id,an_ids{an_ix},erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Exploratory rose plot: Compute and plot ITPC with ERP overlay per condition for single SBJ
    %   Also extracts/plots mean phase angle in T-F window across conditions
%     SBJ05c_ITC_ERP_rose_plot_grp(SBJ_id,conditions,proc_id,an_ids{an_ix},phs_freq_lim,phs_time_lim,erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

