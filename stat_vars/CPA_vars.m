%% Selection Criteria
cpa.diff_id   = 'TarStd';
cpa.elec_list = {'CPz','FCz', 'Fz', 'Cz'};
cpa.time_win  = [0.17 0.4];                 % Epoch to search for differences (in sec)
cpa.min_sig_len = 0.25;                     % Minimum length of significant time (sec)
cpa.ic_rank_max = 11;                       % Maximum index of ICs ranked by relative variance explained


