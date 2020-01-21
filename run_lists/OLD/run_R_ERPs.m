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

%% Run Response-Locked ERPs
proc_id    = 'eeg_full_ft';
an_id      = 'ERP_Z4_R2t1_dm2t0_fl05t20';%,'ERP_Pz_F2t1_dm2t0_fl05t20'
stat_id   = 'DifOut';
plt_id     = 'ts_F2to1_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'off';
fig_ftype  = 'png';
for s = 1:numel(SBJs)
    SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
    SBJ03b_ERP_plot(SBJs{s},stat_id,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    close all;
end
SBJ03c_ERP_plot_grp(SBJs,stat_id,proc_id,an_id,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%     plt_id = 'ts_F2to1_but_evnts_sigPatch';
%     SBJ03c_ERP_plot_grp_butterfly(SBJs,stat_conds{st_ix},proc_id,an_id,plt_id,save_fig,...
%         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
%
%     close all;

