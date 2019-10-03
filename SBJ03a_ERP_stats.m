function SBJ03a_ERP_stats(SBJ,conditions,proc_id,an_id)
% Compute ERPs from preprocessed data:
%   Re-align data to event, select channels and epoch, filter, average, run stats, save
% INPUTS:
%   SBJ [str] - ID of subject to run
%   conditions [str] - label of conditions to compute ERPs for
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
% OUTPUTS:
%   roi_erp [cell] - cell array with outputs of ft_timelockanalysis for each condition
%   stat [ft struct] - output of ft_timelockstatistics

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data 
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Load Data
load([SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat']);

% Select Conditions of Interest
[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(conditions, bhv);

%% Prepare Data
% Realign data to desired event
if ~strcmp(proc.event_type,an.event_type)
    cfg = [];
    % Match desired time to closest sample index
    if strcmp(proc.event_type,'S') && strcmp(an.event_type,'F')
        prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);
        cfg.offset = -(prdm_vars.target + prdm_vars.fb_delay)*clean_trials.fsample;
    elseif strcmp(proc.event_type,'S') && strcmp(an.event_type,'R')
        cfg.offset = -bhv.rt*clean_trials.fsample;
    elseif strcmp(proc.event_type,'F')
        error('F-locked preprocessing can only be used for F-locked analysis!');
    elseif strcmp(proc.event_type,'R')% && strcmp(an.event_type,'S')
        error('Why were you doing R-locked preprocessing?');
        %error('cannot do S-locked analysis with R-locked data!');
    else
        error('unknown combination of proc and an event_types');
    end
    % Convert time axis to new event:
    %   basically: data.time{i} = data.time{i} + offset(i)/data.fsample;
    %   therefore, negative offset will shift time axis "back"
    roi = ft_redefinetrial(cfg, clean_trials);
else
    roi = clean_trials;
end

% Check window consistency
%   Check trial_lim_s is within trial time (round to avoid annoying computer math)
if round(an.trial_lim_s(1)+1/roi.fsample,3) < round(roi.time{1}(1),3) || ...
        round(an.trial_lim_s(2)-1/roi.fsample,3) > round(roi.time{1}(end),3)
    error('an.trial_lim_s is outside data time bounds!');
end
%   Check bsln_lim is within trial_lim_s
if an.bsln_lim(1) < an.trial_lim_s(1) || an.bsln_lim(2) > an.trial_lim_s(2)
    error('an.bsln_lim is outside an.trial_lim_s!');
end

% Select window and channels of interest
cfgs = [];
cfgs.channel = an.ROI;
cfgs.latency = an.trial_lim_s;
roi = ft_selectdata(cfgs, roi);

%% Compute ERPs
% Preprocess the data
cfgpp = [];
cfgpp.hpfilter       = an.hp_yn;
cfgpp.hpfreq         = an.hp_freq;
cfgpp.hpfiltord      = an.hp_filtord; % Leaving blank causes instability error, 1 or 2 works 
cfgpp.lpfilter       = an.lp_yn;
cfgpp.lpfreq         = an.lp_freq;
cfgpp.demean         = an.demean_yn;
cfgpp.baselinewindow = an.bsln_lim;
roi = ft_preprocessing(cfgpp, roi);

% Average ERPs
roi_erp = cell(size(cond_lab));
n_trials = zeros([1 numel(cond_lab)]);
for cond_ix = 1:numel(cond_lab)
    cfg_iavg.trials = find(cond_idx==cond_ix);
    roi_erp{cond_ix} = ft_timelockanalysis(cfg_iavg,roi);
    % Grab n_trials for design matrix
    n_trials(cond_ix) = numel(cfg_iavg.trials);
end

%% Contrast conditions
if strcmp(conditions,'DifOut') || numel(cond_lab)>2
    error(['Analysis not implemented for ' conditions ' yet, too many conditions: ' strjoin(cond_lab,',')]);
%     % Compute Win - Loss
%     cfgdif = [];
%     cfgdif.operation = 'subtract';
%     cfgdif.parameter = 'avg';
%     wn_ix = find(~cellfun(@isempty,strfind(cond_lab,'Wn')));
%     ls_ix = find(~cellfun(@isempty,strfind(cond_lab,'Ls')));
%     difwave = cell(size(wn_ix));
%     for dif_ix = 1:numel(wn_ix)
%         difwave{dif_ix} = roi_erp{wn_ix(dif_ix)};
%         % Compute difference in means (ERP difference wave)
%         difwave{dif_ix}.avg = difwave{dif_ix}.avg - roi_erp{ls_ix(dif_ix)}.avg;
%         difwave{dif_ix} = rmfield(difwave{dif_ix},'trial'); % remove trial since that doesn't make sense anymore
%         % Compute variance of difference in 2 RVs = var(X) + var(Y) - 2*covariance(X,Y)
%         wl_cov = cov(difwave{dif_ix}.avg, roi_erp{ls_ix(dif_ix)}.avg);
%         difwave{dif_ix}.var = difwave{dif_ix}.var + roi_erp{ls_ix(dif_ix)}.var - 2*wl_cov(1,2);
%     end
%     stat_lab = {strrep(cond_lab{wn_ix(1)},'Wn',''), strrep(cond_lab{wn_ix(2)},'Wn','')};
%     
%     cfgdw = [];
%     cfgdw.method = 'stats';
%     cfgdw.statistic = 'ttest2';
%     cfgdw.alpha = 0.05;
%     cfgdw.tail = 0;
%     cfgdw.parameter = 'avg';
%     cfgdw.design = [1 2];
%     stat = ft_timelockstatistics(cfgdw, difwave{:});
end

%% Run Statistics
% Create design matrix
design = zeros(2,sum(n_trials));
for cond_ix = 1:numel(cond_lab)
    if cond_ix==1
        design(1,1:n_trials(cond_ix)) = cond_ix;                                % Conditions (Independent Variable)
        design(2,1:n_trials(cond_ix)) = 1:n_trials(cond_ix);                    % Trial Numbers
    else
        design(1,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= cond_ix; % Conditions (Independent Variable)
        design(2,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= 1:n_trials(cond_ix);
    end
end

% Prepare neighbors layout
% cfgn = [];
% cfgn.method  = 'distance';
% cfgn.layout  = 'ordered';
% cfgn.channel = elecs;
% neighbors    = ft_prepare_neighbours(cfgn,roi_erp_allch{1});
% for ch_ix = 1:numel(roi_erp{1}.label)
%     neighbors(ch_ix).label = roi_erp{1}.label{ch_ix};
%     neighbors(ch_ix).neighblabel = {};
% end

% Calculate statistics
cfg_stat.design = design;
% [stat] = ft_timelockstatistics(cfg_stat, roi_erp{:});
warning('HACK!!!! MEX file work around, just copying roi_erp instead of running real stats!!!');
stat = roi_erp{1};
stat.mask = zeros([numel(stat.label) numel(stat.time)]);

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','roi_erp','stat');

end
