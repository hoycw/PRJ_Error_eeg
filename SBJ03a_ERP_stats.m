function SBJ03a_ERP_stats(SBJ,proc_id,an_id)
% Compute ERPs from preprocessed data:
%   Re-align data to event, select channels and epoch, filter, average, run stats, save
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
% OUTPUTS:
%   roi_erp [cell] - cell array with outputs of ft_timelockanalysis for each condition
%   stat [ft struct] - output of ft_timelockstatistics

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
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
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_final.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

%% Create design matrix
% Select Conditions of Interest
[cond_lab, ~, ~, ~] = fn_condition_label_styles(an.model_lab);
full_cond_idx = fn_condition_index(cond_lab, bhv);
bhv_fields = fieldnames(bhv);
orig_n_trials = numel(bhv.trl_n);
for f_ix = 1:numel(bhv_fields)
    if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
        bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
    end
end

design = cell([1 numel(an.groups)]);
levels = cell([1 numel(an.groups)]);
for grp_ix = 1:numel(an.groups)
    [levels{grp_ix}, ~, ~] = fn_condition_label_styles(an.groups{grp_ix});
    design{grp_ix} = fn_condition_index(levels{grp_ix}, bhv);
end

%% Select Data
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
cfgs.trials  = find(full_cond_idx);
roi = ft_selectdata(cfgs, roi);

%% Prepare Data
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

% Create ANOVA data matrix
cfg = [];
cfg.latency = an.stat_lim;
st_roi = ft_selectdata(cfg, roi);
st_data = nan([numel(bhv.trl_n) numel(st_roi.label) numel(st_roi.time{1})]);
for trl_ix = 1:numel(bhv.trl_n)
    st_data(trl_ix,:,:) = st_roi.trial{trl_ix};
end
if any(isnan(st_data(:))); error('NaN in ANOVA data!'); end

%% Run ANOVA
fprintf('================== Running ANOVA =======================\n');
% Create structure for w2 in fieldtrip style
w2.design    = design;
w2.cond      = an.groups;
w2.time      = st_roi.time{1};
w2.label     = st_roi.label;
w2.dimord    = 'rpt_chan_time';
w2.boot      = zeros([numel(w2.cond) numel(st_roi.label) length(w2.time) an.n_boots]);

% Compute ANOVA and Explained Variance for real model
w2.trial = fn_mass_ANOVA(st_data,design);

% Compute ANOVA for permuted data
rand_design = design;
% b = '';
fprintf('boot #: ');
for boot_ix = 1:an.n_boots
%     m = sprintf(' permutation %d/%d', boot_ix, n_boots);
%     fprintf([b m]); b = repmat('\b',[1 length(m)]);
    fprintf('%i..',boot_ix);
    for grp_ix = 1:numel(an.groups)
        rand_design{grp_ix} = design{grp_ix}(randperm(size(design{grp_ix},1)));
    end
    w2.boot(:,:,:,boot_ix) = fn_mass_ANOVA(st_data,rand_design);
    if mod(boot_ix,20)==0
        fprintf('\n');
    end
end

% Compute statistics
w2.pval = sum(bsxfun(@ge,w2.boot,w2.trial),4)./an.n_boots; % sum(boots>real)/n_boots
w2.zscore   = norminv(1-w2.pval,0,1);
w2.bootmean = mean(w2.boot,4);
w2.bootstd  = std(w2.boot,[],4);
% w2 = rmfield(w2,'boot');
w2.zscore(isinf(w2.zscore)) = norminv(1-1/an.n_boots/2,0,1);

% Multiple Comparisons Correction within Channel
w2.qval = nan(size(w2.pval));
for ch_ix = 1:numel(w2.label)
    [~, ~, ~, w2.qval(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));
end

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','roi','w2','an');

end
