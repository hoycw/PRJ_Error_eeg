% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'ERPEsL';
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.trial_cond  = {'DifFB'};
st.stat_lim    = [-0.025 0.025];            % window in SEC for stats
st.pk_reg_id   = 'sRPE';
st.pk_stat_id  = 'ERPEsL_all_lme_st05t5';            % stat_id from which to get peak time
st.pk_an_id    = 'ERP_Fz_F2t1_dm2t0_fl05t20';            % an_id from which to get peak time
st.measure     = 'mean';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

