function Run_SBJ06(SBJs)
for x = 1: numel(SBJs)
    %ODD02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 'odd_full_ft', 1, 'on', 'ERP_stack_full_events_odd')
    %SBJ02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 1, 'on')
    %SBJ01_preproc(SBJs{x}, 'eeg_full_ft');
    SBJ06a_CPA2(SBJs{x}, 'odd_full_ft', 'ts_F2to13_evnts_sigLine', {'CPz','FCz', 'Fz', 'Cz'}, [0.17, 0.40]);
end
