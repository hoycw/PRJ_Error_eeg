% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'ERPEsL';
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.trial_cond  = {'DifFB'};
st.stat_lim    = [-0.025 0.025];            % window in SEC for stats
st.pk_center   = 0.225;
st.measure     = 'mean';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

