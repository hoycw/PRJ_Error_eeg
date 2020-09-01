% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'SML';
st.z_reg       = 0;
st.trial_cond  = {'DifFB'};
st.stat_lim    = [0.05 0.5];            % window in SEC for stats
st.pk_lim      = [0.1 0.26; 0.18 0.3];
st.pk_sign     = [1; -1];
st.measure     = 'p2p';             % {'p2p', 'mean', 'ts'}
st.grp_method  = 'sbj';             % {'sbj', 'jackknife'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

