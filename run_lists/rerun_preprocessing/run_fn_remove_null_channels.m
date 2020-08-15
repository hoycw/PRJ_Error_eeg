function run_fn_remove_null_channels
    SBJs = {'EEG18', 'EEG19', 'EEG20', 'EEG21', 'EEG22', 'EEG23', 'EEG29'}
    for index = 1:length(SBJs)
        fn_remove_null_channels(SBJs{index}, 'eeg_full_ft')
        fn_remove_null_channels(SBJs{index}, 'odd_full_ft')
    end
end