% Spatial Criteria
cpa.elec_method = {'peak','topo_corr','elec_outlier'};          % consider 'grubbs' as recommended by Wessel?
cpa.topo_SBJ_id = 'goodEEG';
cpa.topo_an_id  = 'ERP_all_S2t1_dm2t0_fl05t20';                   
cpa.topo_lim    = [0.3 0.45];
cpa.topo_cond   = 'Odd';
cpa.topo_pval   = 0.05;
cpa.elec_list = {'CPz','FCz', 'Fz', 'Cz'};
cpa.n_max_elec  = 5;                         % number of peak elecs for method = 'peak'
cpa.min_elec_match = 2;                        % number of electrodes that must match between IC and list
cpa.sd_thresh   = 5;

% Temporal Criteria
cpa.erp_filt  = 0;
cpa.erp_an_id = '';
cpa.diff_id   = 'OddStd';
cpa.time_win  = [0.15 0.4];                 % Epoch to search for differences (in sec)
cpa.min_sig_len = 0.1;                    % Minimum length of significant time (sec)
cpa.alpha       = 0.05;                     % Significance threshold for condition stats

% Variance Criteria
cpa.ic_rank_max = 10;                       % Maximum index of ICs ranked by relative variance explained


