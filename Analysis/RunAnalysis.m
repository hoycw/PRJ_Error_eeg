function runscript(SBJs)
plt_id = 'ts_F15to28_evnts_sigPatch';
odd_plt_id = 'ts_S2to13_evnts_sigPatch';
an_id = 'ERP_Cz_F_trl15t28_flt05t20_stat06';
odd_an_id = 'ERP_Cz_S2t13_flt5t20_st6';
odd_proc_id = 'odd_full_ft';
proc_id = 'eeg_full_ft';
for s = 1:numel(SBJs)
    PLOTSTAT_HitErr_Individ(SBJs{s}, proc_id, plt_id, an_id, 'off', 1, 'png');
    PLOT_4Conds_Individ(SBJs{s}, proc_id, plt_id, an_id, 'off', 1, 'png');
    if startsWith(SBJs{s},'EEG')
        PLOTSTAT_Oddball_Individ(SBJs{s}, odd_proc_id, odd_plt_id, odd_an_id, 'off', 1, 'png');
    end
    %CALC_FRN_Individ(SBJs{s}, proc_id, plt_id, an_id);
end

PLOTSTAT_HitErr_Group(SBJs, plt_id, an_id, 'off', 1, 'png');
PLOT_4Conds_Group(SBJs,  plt_id, an_id, 'off', 1, 'png');
PLOTSTAT_Oddball_Group(SBJs, odd_plt_id, odd_an_id, 'off', 1, 'png');
%PLOTSTATS_FRN_Group(SBJs, proc_id, plt_id, an_id, 'off', 1, '.png');
PLOT_DiffWave_Group(SBJs, plt_id, an_id, 'off', 1, '.png', 'DifWL');
PLOT_DiffWave_Group(SBJs, plt_id, an_id, 'off', 1, '.png', 'DifEH');