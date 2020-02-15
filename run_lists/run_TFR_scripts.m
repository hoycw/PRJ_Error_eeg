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

%% Compute TFRs
proc_id    = 'eeg_full_ft';
%an_ids     = {'TFR_Fz_F2t1_db2t0_fl05t20'};
an_ids     = {'TFR_Fz_F2t1_db2t0_fl1t12b05'};%'TFR_Fz_F2t1_z2t05_fl1t14','TFR_Fz_F2t1_rc2t0_fl1t14','TFR_Fz_F2t1_db2t0_fl1t14'};
erp_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20'};
%an_ids = {'TFR_Fz_F2t1_z2t0_fl1t14','TFR_Pz_F2t1_z2t0_fl1t14'};
conditions = 'DifFB';
plt_id     = 'ts_F2to1_evnts_sigLine';%'ts_F4t1_evnts_sigLine';%
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        SBJ05a_TFR_save(SBJs{s}, proc_id, an_ids{an_ix})
%         SBJ05b_TFR_plot(SBJs{s}, conditions, proc_id, an_ids{an_ix}, plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         SBJ05b_TFR_ERP_plot(SBJs{s}, conditions, proc_id, an_ids{an_ix},erp_ids{an_ix}, plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
%     SBJ05c_TFR_plot_grp(SBJs, conditions, proc_id, an_ids{an_ix}, plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%     SBJ05c_TFR_ERP_plot_grp(SBJs, conditions, proc_id, an_ids{an_ix},erp_ids{an_ix}, plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Compute ITPC
proc_id    = 'eeg_full_ft';
an_ids     = {'ITC_Fz_F2t1_fl05t20'};%{'ITC_Fz_F2t1_fl1t12b05'};%,'ITC_Pz_F2t1_fl1t12b05'};
erp_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20'};
conditions = 'DifFB';
plt_id     = 'ts_F2to1_evnts_sigLine';%'ts_F4t1_evnts_sigLine';%
save_fig    = 1;
fig_vis     = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 2:numel(SBJs)
        SBJ05a_TFR_save(SBJs{s}, proc_id, an_ids{an_ix})
%         SBJ05b_ITC_plot(SBJs{s}, conditions, proc_id, an_ids{an_ix}, plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        SBJ05b_ITC_ERP_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},erp_ids{an_ix},...
            plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
%     SBJ05c_ITC_plot_grp(SBJs, conditions, proc_id, an_ids{an_ix}, plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    SBJ05c_ITC_ERP_plot_grp(SBJs,conditions,proc_id,an_ids{an_ix},erp_ids{an_ix},...
            plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Compute and Plot POW (Time Series)
proc_id    = 'eeg_full_ft';
conditions = 'DifFB';
an_ids     = {'POW_Fz_F2t1_rc2t0_fl4t8','POW_Fz_F2t1_z2t0_fl4t8'};
%an_ids     = {'POW_Fz_F2t1_db2t0_fl4t8','POW_Fz_F2t1_db2t0_fl8t12','POW_Pz_F2t1_db2t0_fl1t4'};
plt_id     = 'ts_F2to1_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        SBJ05a_TFR_save(SBJs{s}, proc_id, an_ids{an_ix});
        
        % Plot evoked POW time series
        SBJ05b_POW_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    SBJ05c_POW_plot_grp(SBJs,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
%     plt_id = 'ts_F2to1_but_evnts_sigPatch';
%     SBJ05c_POW_plot_grp_butterfly(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Compute and Plot POW (Topos)
proc_id    = 'eeg_full_ft';
conditions = 'DifFB';
an_ids     = {'POW_all_F2t1_db2t0_fl4t8','POW_all_F2t1_db2t0_fl8t12','POW_all_F2t1_db2t0_fl1t4'};
topo_plt_ids = {'topo_F18t25','topo_F18t25','topo_F3t45'};
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
%     for s = 1:numel(SBJs)
%         SBJ05a_TFR_save(SBJs{s}, proc_id, an_ids{an_ix});
%         
%         % Plot POW topos
%         SBJ05b_POW_plot_topo_cond(SBJs{s},conditions,proc_id,an_ids{an_ix},...
%             topo_plt_ids{an_ix},save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%         close all;
%     end
    SBJ05c_POW_plot_grp_topo_cond(SBJs,conditions,proc_id,an_ids{an_ix},...
        topo_plt_ids{an_ix},save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end


