%% Time-frequency power and phase analyses for Sequential PE Initial Submission
% Developed over time, but editted 8/20/20 by Colin W Hoy
%   Sup. Fig. 4 and 5: SBJ05c_TFR_ERP_plot_grp

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
SBJ_id = 'goodall';%'good1';
SBJs = fn_load_SBJ_list(SBJ_id);

%% Compute TFRs for Power
proc_id    = 'eeg_full_ft';
an_ids     = {'TFR_Fz_F2t1_db2t0_fl1t12','TFR_Pz_F2t1_db2t0_fl1t12'};
erp_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
conditions = 'DifFB';
plt_id     = 'ts_F2t8_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'svg';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        % Reconstruct and clean raw data, filter, cut trials to event,
        %   baseline correct, select channels, save
        SBJ05a_TFR_save(SBJs{s}, proc_id, an_ids{an_ix});
        
        % Plot TFRs of power data per condition for single SBJ
%         SBJ05b_TFR_plot(SBJs{s}, conditions, proc_id, an_ids{an_ix}, plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot TFRs of power data with ERP overlay per condition for single SBJ
%         SBJ05b_TFR_ERP_plot(SBJs{s}, conditions, proc_id, an_ids{an_ix},erp_ids{an_ix}, plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Plot TFRs of power data per condition for group
%     SBJ05c_TFR_plot_grp(SBJ_id, conditions, proc_id, an_ids{an_ix}, plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Plot TFRs of power data with ERP overlay per condition for group
    %*** Sup. Fig. 4 (Fz) and 5 (Pz)
    SBJ05c_TFR_ERP_plot_grp(SBJ_id, conditions, proc_id, an_ids{an_ix},erp_ids{an_ix}, plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Compute Phase for ITPC
proc_id    = 'eeg_full_ft';
an_ids     = {'PHS_Fz_F2t1_fl1t12','PHS_Pz_F2t1_fl1t12'};
erp_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
conditions = 'DifFB';
plt_id     = 'ts_F2to1_evnts_sigLine';
save_fig    = 1;
fig_vis     = 'on';
fig_ftype  = 'png';

% FRN Phase ROI
phs_freq_lim = [1 4];
phs_time_lim = [0.3 0.33];

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        % Reconstruct and clean raw data, filter, cut trials to event,
        %   baseline correct, select channels, save
        SBJ05a_TFR_save(SBJs{s}, proc_id, an_ids{an_ix});
        
        % Compute and plot ITPC with ERP overlay per condition for single SBJ
%         SBJ05b_ITC_ERP_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Exploratory rose plot: Compute and plot ITPC with ERP overlay per condition for single SBJ
        %   Also extracts/plots mean phase angle in T-F window across conditions
%         SBJ05b_ITC_ERP_rose_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},phs_freq_lim,phs_time_lim,erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Compute and plot ITPC with ERP overlay per condition for group
%     SBJ05c_ITC_ERP_plot_grp(SBJ_id,conditions,proc_id,an_ids{an_ix},erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % Exploratory rose plot: Compute and plot ITPC with ERP overlay per condition for single SBJ
    %   Also extracts/plots mean phase angle in T-F window across conditions
%     SBJ05c_ITC_ERP_rose_plot_grp(SBJ_id,conditions,proc_id,an_ids{an_ix},phs_freq_lim,phs_time_lim,erp_ids{an_ix},...
%             plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

