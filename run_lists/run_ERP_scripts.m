%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%%
addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Run TT ERPs
SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP19',...
           'EEG01','EEG02','EEG03','EEG04','EEG06','EEG08'};

conditions = 'Out';
proc_id    = 'eeg_full_ft';
an_id      = 'ERP_Cz_F25to1_flt05to20_st06';
plt_id     = 'ts_F25to1_evnts_sigPatch';
save_fig   = 0;
fig_vis    = 'on';
fig_ftype  = 'png';
for s = 1:numel(SBJs)
%     SBJ03a_ERP_stats(SBJs{s},conditions,proc_id,an_id);
    SBJ03b_ERP_plot_stats(SBJs{s},conditions,proc_id,an_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end

% 
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
odd_SBJs = {'EEG01','EEG02','EEG03','EEG04','EEG06','EEG08'};

proc_id   = 'odd_full_ft';
an_id     = 'ERP_Cz_S2t13_flt5t20_st6';
plt_id    = 'ts_S2to13_evnts_sigPatch';
fig_vis   = 'on';
save_fig  = 1;
fig_ftype = 'png';

for s = 4:numel(SBJs)
    oddball_stats_plot(SBJs{s}, proc_id, plt_id, an_id, fig_vis, save_fig, fig_ftype)
end