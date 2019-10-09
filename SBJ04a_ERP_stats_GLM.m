function SBJ04a_ERP_stats_GLM(SBJ,proc_id,an_id,stat_id)
error('load roi dont compute');
% Compute ERPs from preprocessed data:
%   Re-align data to event, select channels and epoch, filter, average, run stats, save
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   roi [ft struct] - preprocessed data (can be averaged to get ERP)
%   beta [pseudo-ft struct] - output of GLM

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
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Load Data
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_final.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

%% Create Design Matrix
% Select Conditions of Interest
[cond_lab, ~, ~, ~] = fn_condition_label_styles(st.model_lab);
full_cond_idx = fn_condition_index(cond_lab, bhv);
bhv_fields = fieldnames(bhv);
orig_n_trials = numel(bhv.trl_n);
for f_ix = 1:numel(bhv_fields)
    if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
        bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
    end
end

design = nan([numel(bhv.trl_n) numel(st.regressors)]);
for r_ix = 1:numel(st.regressors)
    design(:,r_ix) = fn_build_regressor(st.regressors{r_ix}, bhv);
end

% levels = cell([1 numel(st.groups)]);
% for grp_ix = 1:numel(st.groups)
%     [levels{grp_ix}, ~, ~] = fn_condition_label_styles(st.groups{grp_ix});
%     design{grp_ix} = fn_condition_index(levels{grp_ix}, bhv);
% end

%% Select Data for ERP
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

% Select window and channels of interest
cfgs = [];
cfgs.channel = an.ROI;
cfgs.latency = an.trial_lim_s;
cfgs.trials  = find(full_cond_idx);
roi = ft_selectdata(cfgs, roi);

%% Preprocess Data for ERP
cfgpp = [];
cfgpp.hpfilter       = an.hp_yn;
cfgpp.hpfreq         = an.hp_freq;
cfgpp.hpfiltord      = an.hp_filtord; % Leaving blank causes instability error, 1 or 2 works 
cfgpp.lpfilter       = an.lp_yn;
cfgpp.lpfreq         = an.lp_freq;
cfgpp.demean         = an.demean_yn;
cfgpp.baselinewindow = an.bsln_lim;
roi = ft_preprocessing(cfgpp, roi);

%% Create ANOVA data matrix
cfg = [];
cfg.latency = st.stat_lim;
st_roi = ft_selectdata(cfg, roi);
st_data = nan([numel(bhv.trl_n) numel(st_roi.label) numel(st_roi.time{1})]);
for trl_ix = 1:numel(bhv.trl_n)
    st_data(trl_ix,:,:) = st_roi.trial{trl_ix};
end
if any(isnan(st_data(:))); error('NaN in ANOVA data!'); end

%% Run ANOVA
fprintf('================== Running GLM =======================\n');
% Create structure for beta in fieldtrip style
beta.design = design;
beta.cond   = st.regressors;
beta.time   = st_roi.time{1};
beta.label  = st_roi.label;
beta.dimord = 'rpt_chan_time';
beta.boot   = zeros([numel(beta.cond) numel(st_roi.label) length(beta.time) st.n_boots]);

% Compute ANOVA and Explained Variance for real model
beta.trial = fn_mass_GLM(design,st_data,0);

% Generate randomization vectors
rand_vecs = zeros([size(design,1) numel(st.n_boots)]);
fprintf('Generating randomization vectors: ');
for boot_ix = 1:st.n_boots
    if mod(boot_ix,50)==0
        fprintf('%i..',boot_ix);
    end
    rand_vecs(:,boot_ix) = randperm(size(design,1));
end
fprintf('\n');

% Compute ANOVA for permuted data (permute one regressor at a time)
rand_design = design;
for grp_ix = 1:numel(st.regressors)
    if ~strcmp(st.regressors{grp_ix},'off')
        fprintf('%s (%i / %i) boot #: ',st.regressors{grp_ix},grp_ix,numel(st.regressors));
        for boot_ix = 1:st.n_boots
            fprintf('%i..',boot_ix);
            rand_design(:,grp_ix) = design(rand_vecs(:,boot_ix),grp_ix);
            beta.boot(:,:,:,boot_ix) = fn_mass_GLM(rand_design,st_data,0);
            if mod(boot_ix,20)==0
                fprintf('\n');
            end
        end
    end
end

% Compute statistics
beta.pval = sum(bsxfun(@ge,beta.boot,beta.trial),4)./st.n_boots; % sum(boots>real)/n_boots
beta.zscore   = norminv(1-beta.pval,0,1);
beta.bootmean = mean(beta.boot,4);
beta.bootstd  = std(beta.boot,[],4);
% beta = rmfield(beta,'boot');
beta.zscore(isinf(beta.zscore)) = norminv(1-1/st.n_boots/2,0,1);

% Multiple Comparisons Correction within Channel
beta.qval = nan(size(beta.pval));
for ch_ix = 1:numel(beta.label)
    for r_ix = 1:numel(st.regressors)
        [~, ~, ~, beta.qval(r_ix,ch_ix,:)] = fdr_bh(squeeze(beta.pval(r_ix,ch_ix,:)));
    end
end

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',stat_id,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','bhv','roi','beta');

end
