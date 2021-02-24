% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'ERPEsL';
st.model_cond  = 'DifFB';           % trial conditions to model
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.stat_cond   = 'DifFB';
st.stat_lim    = [0.05 0.5];            % window in SEC for stats
st.measure     = 'ts';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

% Add bias to win probability to approximate subjective ratings
st.bias_reg    = {'pWin'};
st.bias        = [0.25];
st.bias_cond   = {'Hd'};

st.model_id    = [st.model_lab '_pW25hd_' st.model_cond];
warning(['Manually setting bias in model_id: ' st.model_id]);

