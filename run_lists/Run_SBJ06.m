function Run_SBJ06(SBJs)
for x = 1: numel(SBJs)
    SBJ06a_CPA2(SBJs{x}, 'odd_full_ft', 'ts_F2to13_evnts_sigLine', {'FCz', 'Fz', 'Cz'}, [0.17, 0.40]);
    %SBJ06a_CPA2(SBJs{x}, 'odd_full_ft', 'ts_F2to13_evnts_sigLine', {'CPz', 'Pz', 'Oz'}, [0.25, 0.40]);
    %SBJ06a_CPA_TT(SBJs{x}, 'eeg_full_ft', 'ts_F2to13_evnts_sigLine', {'FCz', 'Cz', 'Fz'}, [0.17, 0.40]);
    %SBJ06a_CPA_TT(SBJs{x}, 'eeg_full_ft', 'ts_F2to13_evnts_sigLine', {'CPz', 'Pz', 'Oz'}, [0.25, 0.40]);
end
