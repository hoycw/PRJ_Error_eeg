% Stat Parameters
st.an_style    = 'lme';
st.model_lab   = 'sRL';
st.z_reg       = 1;                 % 0/1: zscore regressors?
st.trial_cond  = {'DifFB'};
st.stat_lim    = [-0.05 0.05];            % window in SEC for stats
st.pk_trial_cond = 'All';                       % Which trial_cond set of ERPs to load
st.pk_erp_cond = 'All';                         % Which condition to take ERP peak from
st.pk_lim      = [0.18 0.3];                    % Window to search for peak
st.pk_sign     = -1;                            % Sign of peak to find
st.pk_an_id    = 'ERP_Fz_F2t1_dm2t0_fl05t20';            % an_id from which to get peak time
st.measure     = 'mean';             % {'ts', 'p2p', 'mean'}
st.n_boots     = 1000;             % Repetitions for non-parametric stats
st.alpha       = 0.05;
st.mcp_method  = 'FDR';

% % Fieldtrip stats
% cfg_stat = [];
% cfg_stat.latency          = an.stat_lim;
% cfg_stat.channel          = 'all';
% cfg_stat.parameter        = 'trial';
% cfg_stat.method           = 'montecarlo';
% cfg_stat.statistic        = 'ft_statfun_indepsamplesT';
% cfg_stat.correctm         = 'cluster';
% cfg_stat.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
% cfg_stat.clusterstatistic = 'maxsum';
% cfg_stat.clustertail      = 0;
% cfg_stat.tail             = 0; %two sided
% cfg_stat.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
% cfg_stat.alpha            = 0.05;
% cfg_stat.numrandomization = an.n_boots;
% cfg_stat.neighbours       = [];%neighbors;
% % cfg_stat.minnbchan        = 0;
% cfg_stat.ivar             = 1;  %row of design matrix containing independent variable
% % cfg_stat.uvar             = 2;  %row containing dependent variable, not needed for indepsamp

% % Fieldtrip averaging within SBJ
% cfg_iavg = [];
% cfg_iavg.keeptrials = 'no';
% 
% % Fieldtrip averaging across SBJs
% cfg_gavg = [];
% cfg_gavg.channel        = {'all'};
% cfg_gavg.latency        = an.stat_lim;
% cfg_gavg.keepindividual = 'no';
% cfg_gavg.method         = 'across';
% cfg_gavg.parameter      = 'avg';
% 
% % cfg_tfr.method = 'wavelet'
% % cfg_tfr.output = 'pow';
% % cfg_tfr.taper = 'hanning';
% % cfg_tfr.foi = [2:30]; 
% % cfg_tfr.width = 2; %default
% % cfg_tfr.toi  = 0:0.004:2.8; 
