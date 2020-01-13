% Spatial Criteria
st.elec_list = {'CPz','FCz', 'Fz', 'Cz'};
st.elec_method = 'peak';                    % consider 'grubbs' as recommended by Wessel?
st.n_max_elec  = 5;                         % number of peak elecs for method = 'peak'
st.min_elec_match = 2;                        % number of electrodes that must match between IC and list

% Temporal Criteria
st.diff_id   = 'TarStd';
st.time_win  = [0.17 0.4];                 % Epoch to search for differences (in sec)
st.min_sig_len = 0.01;                    % Minimum length of significant time (sec)
st.alpha       = 0.05;                     % Significance threshold for condition stats

% Variance Criteria
st.ic_rank_max = 11;                       % Maximum index of ICs ranked by relative variance explained


