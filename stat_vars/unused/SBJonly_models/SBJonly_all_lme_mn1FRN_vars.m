% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'SBJonly';
st.model_cond  = 'DifFB';           % trial conditions to model
st.model_id    = [st.model_lab '_' st.model_cond];
st.z_reg       = 0;
st.stat_cond   = {'DifFB'};
st.stat_lim    = [-0.05 0.05];            % window in SEC for stats
st.pk_stat_cond = 'All';                       % Which stat_cond set of ERPs to load
st.pk_erp_cond = 'All';                         % Which condition to take ERP peak from
st.pk_lim      = [0.18 0.3];                    % Window to search for peak
st.pk_sign     = -1;                            % Sign of peak to find
st.pk_an_id    = 'ERP_Fz_F2t1_dm2t0_fl05t20';   % an_id from which to get peak
st.measure     = 'mean';             % {'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

