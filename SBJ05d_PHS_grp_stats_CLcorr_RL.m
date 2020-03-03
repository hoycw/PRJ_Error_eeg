function SBJ05d_PHS_grp_stats_CLcorr_RL(SBJs,proc_id,an_id,stat_id)
% Run Circular-Linear correlation for phase at each time-frequency point
% across all SBJ (no random effects...), then Fisher z-transformed r values
% and t-test against 0 (and FDR correct)
%   Only for one channel now...
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults
addpath([app_dir 'CircStat/']);

%% Load Data 
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~an.itpc; error('why run this without ITPC an_vars?'); end

stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.an_style,'CLcorr'); error('stat_id not using circular-linear correlation!'); end

% Select conditions (and trials)
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Behavior
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
for s = 1:numel(SBJs)
    % Load data
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    bhvs{s} = tmp.bhv;
    
    % Select Conditions of Interest
    full_cond_idx{s} = fn_condition_index(cond_lab, bhvs{s});
    bhv_fields = fieldnames(bhvs{s});
    orig_n_trials = numel(bhvs{s}.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhvs{s}.(bhv_fields{f_ix}))==orig_n_trials
            bhvs{s}.(bhv_fields{f_ix}) = bhvs{s}.(bhv_fields{f_ix})(full_cond_idx{s}~=0);
        end
    end
    n_trials(s) = numel(bhvs{s}.trl_n);
    
    clear tmp
end

%% Load Data and Build Model
cfgs  = []; cfgs.latency = st.stat_lim;
model = zeros([sum(n_trials) numel(reg_lab)]);
sbj_factor  = zeros([sum(n_trials) 1]);
for s = 1:numel(SBJs)
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    % Load behavior
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_' proc_id '_' an_id '.mat'],'tfr');
    if numel(tfr.label)>1; error('assuming single channel for now!'); end
    
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_tfr = ft_selectdata(cfgs, tfr);
    
    
    % Initialize matrices now that we know time-frequency axes
    if s==1
        time_vec  = st_tfr.time;
        fois      = st_tfr.freq;
        if strcmp(st.measure,'ts')
            data  = nan([sum(n_trials) numel(fois) numel(time_vec)]);
        elseif strcmp(st.measure,'mean')
            data  = nan([sum(n_trials) numel(fois)]);
        else; error(['unknown st.measure: ' st.measure]);
        end
    end
       
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = st_tfr.time;
        fois     = st_tfr.freq;
        if strcmp(st.measure,'ts')
            data  = nan([sum(n_trials) numel(fois) numel(time_vec)]);
        elseif strcmp(st.measure,'mean')
            data  = nan([sum(n_trials) numel(fois)]);
        else; error(['unknown st.measure: ' st.measure]);
        end
        
        % Load RL Model
        tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
        % model(1:n_trials(s),:) = tmp.model;
        % Z-score SBJ model regressors
        sbj_model = NaN(size(tmp.model));
        if st.z_reg
            for reg_ix = 1:numel(reg_lab)
                sbj_model(:,reg_ix) = ...
                    (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
            end
        else
            sbj_model = model;
        end
        model(1:n_trials(s),:) = sbj_model;
        
        % Load and add data
squeeze(angle(st_tfr.fourierspctrm(:,1,f_ix,t_ix)))
        if strcmp(st.measure,'ts')
            data(1:n_trials(s),:,:) = squeeze(st_tfr.powspctrm);
        elseif strcmp(st.measure,'mean')
            data(1:n_trials(s),:) = squeeze(mean(st_tfr.powspctrm,4));
        else; error(['unknown st.measure: ' st.measure]);
        end
        
        % Track SBJ
        sbj_factor(1:n_trials(s),end) = s*ones([n_trials(s) 1]);
    else
        % Load RL Model
        tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
        % model(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),:) = tmp.model;
        % Z-score SBJ model regressors
        sbj_model = NaN(size(tmp.model));
        if st.z_reg
            for reg_ix = 1:numel(reg_lab)
                sbj_model(:,reg_ix) = ...
                    (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
            end
        else
            sbj_model = model;
        end
        model(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),:) = sbj_model;
        
        % Load and add data
        if strcmp(st.measure,'ts')
            data(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),:,:) = squeeze(st_tfr.powspctrm);
        elseif strcmp(st.measure,'mean')
            data(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),:) = squeeze(mean(st_tfr.powspctrm,4));
        else; error(['unknown st.measure: ' st.measure]);
        end
        
        % Track SBJ
        sbj_factor(sum(n_trials(1:s-1))+1:sum(n_trials(1:s))) = s*ones([n_trials(s) 1]);
    end
    
    clear bhv tmp tfr st_tfr model full_cond_idx trial_idx_perm
end
fprintf('\t\t SBJ permutation stats complete:');
toc

%% Build table
fprintf('========================== Running Stats ==========================\n');
% Stats for each regressor and TFR point
        phs_corr  = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec)]);
        phs_null  = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec) st.n_boots]);
    % Phase-Model Correlations
    trial_idx_perm = zeros([size(model,1) st.n_boots]);
    for reg_ix = 1:numel(reg_lab)
        fprintf('%s (%d/%d) freq: ',reg_lab{reg_ix},reg_ix,numel(reg_lab));
        for f_ix = 1:numel(fois)
            fprintf('%.3f..',fois(f_ix));
            for t_ix = 1:numel(time_vec)
                [phs_corr(s,reg_ix,f_ix,t_ix), ~] = circ_corrcl(...
                    data, model(:,reg_ix));
                
                % Generate null distribution
                for boot_ix = 1:st.n_boots
                    if reg_ix==1 && f_ix==1 && t_ix==1
                        trial_idx_perm(:,boot_ix) = randperm(size(model,1));
                    end
                    [phs_null(s,reg_ix,f_ix,t_ix,boot_ix), ~] = circ_corrcl(...
                        squeeze(angle(st_tfr.fourierspctrm(:,1,f_ix,t_ix))), ...
                        model(trial_idx_perm(:,boot_ix),reg_ix));
                end
                % Z-score correlation values
                phs_zcorr(s,reg_ix,f_ix,t_ix) = ...
                    (phs_corr(s,reg_ix,f_ix,t_ix)-mean(phs_null(s,reg_ix,f_ix,t_ix,:),5)) / ...
                    std(phs_null(s,reg_ix,f_ix,t_ix,:),[],5);
            end
        end
        fprintf('\n');
    end


pvals = nan([numel(reg_lab) numel(fois) numel(time_vec)]);
for reg_ix = 1:numel(reg_lab)
    for f_ix = 1:numel(fois)
        for t_ix = 1:numel(time_vec)
            [~, pvals(reg_ix,f_ix,t_ix)] = ttest(squeeze(phs_zcorr(:,reg_ix,f_ix,t_ix)));
            % Previously: ttest(atanh(squeeze(phs_corr(:,reg_ix,f_ix,t_ix))));
        end
    end
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2)*size(pvals,3) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2) size(pvals,3)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

%% Save Results
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(stat_out_dir,'dir')
    [~] = mkdir(stat_out_dir);
end
stat_out_fname = [stat_out_dir 'GRP_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','phs_corr','phs_zcorr','qvals','SBJs');%phs_null

end
