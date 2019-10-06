% Time Parameters
st.ep_lab   = 'S';
st.evnt_lab = 'S';
st.stat_lim = [0 1.5];
st.lim_adj  = {'', ''};
st.cust_win = 0;            % custom windows per trials
st.min_rt   = 0.35;
st.alpha    = 0.05;

% Sliding Window Parameters (in sec)
st.win_len  = 0.1;
st.win_step = 0.025;

% ANOVA Parameters
st.model_lab   = 'CNI_p_PC';
st.regress_rt  = 0;    % regrees reaction time off before running ANOVA
st.groups      = {'CNI', 'pCNI', 'PC'};
st.trial_cond  = {'all'};
st.n_boots     = 1000;

% RT Correlation Parameters
st.rt_corr       = 0;
cfg_rt = [];
cfg_rt.parameter        = 'powspctrm';
cfg_rt.statistic        = 'ft_statfun_correlationT';
cfg_rt.method           = 'montecarlo';
cfg_rt.numrandomization = st.n_boots;
cfg_rt.correctm         = 'cluster';
cfg_rt.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_rt.clusterstatistic = 'maxsum';
cfg_rt.clustertail      = 0;
cfg_rt.tail             = 0; %two sided
cfg_rt.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_rt.computestat      = 'yes';
cfg_rt.computeprob      = 'yes';
cfg_rt.alpha            = st.alpha;
cfg_rt.neighbours       = [];

