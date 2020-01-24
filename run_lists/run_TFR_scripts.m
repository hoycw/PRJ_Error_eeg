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
           'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG10','EEG12'};
% Bad SBJ:
%   EP01, EP02, EP05- recording errors
%   EP03- low quality
%   EP09- ???
%   EP12- ???
%   EP13- ???
%   EEG05- ???
%   EEG09- multiple blocks, needs redo???
%   EEG11- recording failure

%% Compute TFRs
proc_id    = 'eeg_full_ft';
an_ids     = {'TFR_Fz_F2t1_z2t0_fl1t14','TFR_Fz_F2t1_z2t0_fl1t14'};
conditions = 'DifFB';
% plt_id    = 'ts_F4t1_evnts_sigLine';%'ts_F2to1_evnts_sigLine';
% save_fig    = 1;
% fig_vis     = 'off';
% fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
        SBJ05a_TFR_save(SBJs{s}, proc_id, an_ids{an_ix})
    end
end
%      SBJ05b_TFR_plot(SBJs{s}, 'DifOut', proc_id, an_id, save_fig);
%      SBJ05b_TFR_plot(SBJs{s}, 'Out', proc_id, an_id, save_fig);

%% Compute TFRs
for s = 1:numel(SBJs)
    SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
    SBJ03b_ERP_plot(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

%% Plot POW Topos
proc_id    = 'eeg_full_ft';
conditions = 'DifFB';
an_id      = 'POW_all_F2t1_dm2t0_fl4t8';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

% for s = 1:numel(SBJs)
% %     SBJ05a_TFR_save(SBJs{s}, proc_id, an_id);
%     SBJ05b_TFR_plot(SBJs{s}, conditions, proc_id, an_ids, save_fig);
% end
% 
% SBJ05c_TFR_plot_grp(SBJs,conditions,proc_id,an_ids,save_fig);
% 
% for s = 1:numel(SBJs)
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
% end
% 
% % FRN Group plot
% plt_id    = 'topo_F18t25';
% SBJ03c_ERP_plot_grp_topo_cond(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% 
% % P3 Group Plot
% plt_id    = 'topo_F3t45';
% SBJ03c_ERP_plot_grp_topo_cond(SBJs,conditions,proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
% 
% 
% 
% % Full time series by condition
% 
% 
