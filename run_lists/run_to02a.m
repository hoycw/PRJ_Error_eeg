function run_to02a(SBJs)
for ix = 1:numel(SBJs)
    SBJ01_preproc(SBJs{ix}, 'eeg_full_ft')
    %SBJ02a_artifact_rejection(SBJs{ix}, 'eeg_full_ft', 1, 'on')
    %SBJ01_preproc(SBJs{ix}, 'odd_full_ft')
    %ODD02a_artifact_rejection(SBJs{ix}, 'eeg_full_ft', 'odd_full_ft', 1, 'on', 'ERP_stack_full_events_odd')
end
