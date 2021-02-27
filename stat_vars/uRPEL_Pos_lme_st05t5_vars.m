% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'uRPEL';
st.model_cond  = 'DifFB';           % trial conditions to model
st.model_id    = [st.model_lab '_' st.model_cond];
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.stat_cond   = 'Pos';
st.stat_lim    = [0.05 0.5];            % window in SEC for stats
st.measure     = 'ts';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

