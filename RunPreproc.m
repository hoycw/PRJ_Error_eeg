function RunPreproc(SBJs)
proc_id = 'eeg_full_ft';
odd_proc_id = 'odd_full_ft';
gen_figs = 1;
fig_vis = 1;
for ix = 1:numel(SBJs)
    SBJ02a_artifact_rejection(SBJs{ix}, proc_id, gen_figs, 'on');
end