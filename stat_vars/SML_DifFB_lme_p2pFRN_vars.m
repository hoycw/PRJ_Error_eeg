% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'SML';
st.model_cond  = 'DifFB';           % trial conditions to model
st.model_id    = [st.model_lab '_' st.model_cond];
st.z_reg       = 1;
st.stat_cond   = 'DifFB';
st.stat_lim    = [0.05 0.5];            % window in SEC for stats
st.measure     = 'p2p';             % {'p2p', 'mean', 'ts'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

% Peak time selection parameters
st.pk_lim      = [0.1 0.26; 0.18 0.3];
st.pk_sign     = [1; -1];
st.grp_method  = 'sbj';             % {'sbj', 'jackknife'}

