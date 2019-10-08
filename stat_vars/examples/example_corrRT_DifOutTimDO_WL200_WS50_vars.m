% Model Factors and Levels
model_lab   = 'corrRT_DifOutTimDO';

% ANOVA Parameters
regress_rt  = 0;    % regrees reaction time off before running ANOVA
groups      = {'Dif', 'Out', 'Tim'};
levels      = {{'Ez','Hd'},{'Wn','Ls'},{'Er','Lt'}};
anova_terms = [1 0 0; 0 1 0; 0 0 1; 1 1 0];

% Stats parameters
n_boots = 1000;

% Sliding Window Parameters
win_len = 200;
win_step = 50;

% RT Correlation Parameters
rt_correlation = 1;
cfg_rt = [];
cfg_rt.parameter        = 'powspctrm';
cfg_rt.statistic        = 'ft_statfun_correlationT';
cfg_rt.method           = 'montecarlo';
cfg_rt.numrandomization = n_boots;
cfg_rt.correctm         = 'cluster';
cfg_rt.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_rt.clusterstatistic = 'maxsum';
cfg_rt.clustertail      = 0;
cfg_rt.tail             = 0; %two sided
cfg_rt.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_rt.computestat      = 'yes';
cfg_rt.computeprob      = 'yes';
cfg_rt.alpha            = 0.05;
cfg_rt.neighbours       = [];

