% Data selection
an.ROI         = {'Cz'};             % Channel to be analyzed
an.event_type  = 'F';           % event around which to cut trials
an.trial_lim_s = [-0.2 1];       % window in SEC for cutting trials

% ERP Filtering
an.demean_yn   = 'yes';
an.bsln_lim    = [-0.2 0];    % window in SEC for baseline correction
an.lp_yn       = 'yes';
an.lp_freq     = 20;
an.hp_yn       = 'yes';
an.hp_freq     = 0.5;
an.hp_filtord  = 4;

% ANOVA Parameters
an.model_lab   = 'Out';
an.groups      = {'Out'};
an.trial_cond  = {'all'};
an.stat_lim    = [0 0.6];            % window in SEC for stats
an.n_boots     = 1000;             % Repetitions for non-parametric stats
an.alpha       = 0.05;
% an.anova_terms = [1 0 0; 0 1 0; 0 0 1; 1 1 0];

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
