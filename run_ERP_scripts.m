%% Run Oddball ERPs
SBJs = {'EEG01','EEG02','EEG03','EEG04','EEG06','EEG08'};

proc_id   = 'odd_full_ft';
an_id     = 'ERP_Cz_S2t13_flt5t20_st6';
plt_id    = 'ts_S2to13_evnts_sigPatch';
fig_vis   = 'on';
save_fig  = 1;
fig_ftype = 'png';

for s = 4:numel(SBJs)
    oddball_stats_plot(SBJs{s}, proc_id, plt_id, an_id, fig_vis, save_fig, fig_ftype)
end