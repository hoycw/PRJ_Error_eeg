% Spatial Criteria
cpa.elec_list = {'CPz','FCz', 'Fz', 'Cz'};
cpa.elec_method = 'peak';                    % consider 'grubbs' as recommended by Wessel?
cpa.n_max_elec  = 5;                         % number of peak elecs for method = 'peak'
cpa.min_elec_match = 2;                        % number of electrodes that must match between IC and list

% Temporal Criteria
cpa.diff_id   = 'TarStd';
cpa.time_win  = [0.17 0.4];                 % Epoch to search for differences (in sec)
cpa.min_sig_len = 0.01;                    % Minimum length of significant time (sec)
cpa.alpha       = 0.05;                     % Significance threshold for condition stats

% Variance Criteria
cpa.ic_rank_max = 11;                       % Maximum index of ICs ranked by relative variance explained


