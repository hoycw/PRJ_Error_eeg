% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'EsRPEL';
st.model_cond  = 'DifFB';           % trial conditions to model
st.model_id    = [st.model_lab '_' st.model_cond];
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.stat_cond   = 'DifFB';
st.stat_lim    = [-0.025 0.025];            % window in SEC for stats
st.measure     = 'erp_mean';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

% Peak time selection parameters
st.pk_reg_id   = 'Lik';
st.pk_stat_id  = 'ERPEsL_DifFB_lme_st05t5';            % stat_id from which to get peak time
st.pk_an_id    = 'ERP_Fz_F2t1_dm2t0_fl05t20';            % an_id from which to get peak time

