function SBJ05d_PHS_grp_stats_zCLcorr_RL(SBJs,proc_id,an_id,stat_id)
% Run Circular-Linear correlation for phase at each time-frequency point
% within SBJ, z-score that value to null distribution, then t-test those
% z-scored correlation values vs. 0 across SBJ (and FDR correct)
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
if ~an.complex; error('why run this without ITPC an_vars?'); end

stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.an_style,'zCLcorr'); error('stat_id not using z-scored circular-linear correlation!'); end

% Select conditions (and trials)
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Data and Build Model
cfgs  = []; cfgs.latency = st.stat_lim;
tic
for s = 1:numel(SBJs)
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    % Load behavior
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
    % model(1:n_trials(s),:) = tmp.model;
    % Z-score SBJ model regressors
    model = NaN(size(tmp.model));
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            model(:,reg_ix) = ...
                (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
        end
    else
        model = tmp.model;
    end
    
    % Load data
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_' proc_id '_' an_id '.mat'],'tfr');
    if numel(tfr.label)>1; error('assuming single channel for now!'); end
    
    % Select time and trials of interest
    model = model(full_cond_idx,:);
    cfgs.trials  = find(full_cond_idx);
    st_tfr = ft_selectdata(cfgs, tfr);
    
    % Initialize matrices now that we know time-frequency axes
    if s==1
        time_vec  = st_tfr.time;
        fois      = st_tfr.freq;
        phs_corr  = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec)]);
        phs_zcorr = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec)]);
        phs_null  = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec) st.n_boots]);
    end
    
    % Phase-Model Correlations
    trial_idx_perm = zeros([size(model,1) st.n_boots]);
    for reg_ix = 1:numel(reg_lab)
        fprintf('%s (%d/%d) freq: ',reg_lab{reg_ix},reg_ix,numel(reg_lab));
        for f_ix = 1:numel(fois)
            fprintf('%.3f..',fois(f_ix));
            for t_ix = 1:numel(time_vec)
                [phs_corr(s,reg_ix,f_ix,t_ix), ~] = circ_corrcl(...
                    squeeze(angle(st_tfr.fourierspctrm(:,1,f_ix,t_ix))), model(:,reg_ix));
                
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
    fprintf('End of SBJ %d / %d: %s\n',s,numel(SBJs),SBJs{s});
    toc
    
%     % potential parallelized:
%     % Extract phase angle and reshape matrix so it is nTrials*nFeatures
%     phs = reshape(angle(st_tfr.fourierspctrm(:,1,:,:)),[phs_dims(1) prod(phs_dims(2:end))]);
%     [rho, ~] = circ_corrcl(phs, model(:,reg_ix));
%     final_rho = reshape(rho,dims(2:end));
    
    clear bhv tmp tfr st_tfr model full_cond_idx trial_idx_perm
end
fprintf('\t\t SBJ permutation stats complete:');
toc

%% Build table
fprintf('========================== Running Stats ==========================\n');
% Stats for each regressor and TFR point
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
