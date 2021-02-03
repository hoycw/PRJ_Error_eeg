% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'VML';
st.model_cond  = 'DifFB';           % trial conditions to model
st.model_id    = [st.model_lab '_' st.model_cond];
st.z_reg       = 0;
st.stat_cond   = 'DifFB';
st.stat_lim    = [0.05 0.5];            % window in SEC for stats
st.measure     = 'ts';             % {'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

