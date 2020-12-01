% Stat Parameters
st.an_style    = 'reg';
st.model_lab   = 'P3aP3b_grpMW1';
st.model_cond  = 'All';           % trial conditions to model
st.model_id    = [st.model_lab '_' st.model_cond];
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.stat_cond   = 'All';
st.stat_lim    = [-0.05 0.05];            % window in SEC for stats
st.measure     = 'erp_mean';             % {'ts', 'p2p', 'mean', 'erp_mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

% Peak time selection parameters
st.pk_reg_id   = 'uRPE';
st.pk_stat_id  = 'ERPEsL_all_lme_st05t5';            % stat_id from which to get peak time
st.pk_an_id    = 'ERP_Pz_F2t1_dm2t0_fl05t20';            % an_id from which to get peak time

