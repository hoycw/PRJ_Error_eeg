% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'ERPEsL';
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.trial_cond  = {'DifFB'};
st.stat_lim    = [-0.05 0.05];            % window in SEC for stats
st.pk_trial_cond = 'All';                       % Which trial_cond set of ERPs to load
st.pk_erp_cond = 'All';                         % Which condition to take ERP peak from
st.pk_lim      = [0.25 0.5];                    % Window to search for peak
st.pk_sign     = 1;                            % Sign of peak to find
st.pk_an_id    = 'ERP_Pz_F2t1_dm2t0_fl05t20';            % an_id from which to get peak time
st.measure     = 'mean';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

