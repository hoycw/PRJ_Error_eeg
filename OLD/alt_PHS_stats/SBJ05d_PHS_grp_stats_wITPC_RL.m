function SBJ05d_PHS_grp_stats_wITPC_RL(SBJ_id,proc_id,an_id,stat_id)
error('Do NOT use this analysis, use SBJ05d_PHS_grp_stats_ITPC_jkLME_RL instead!');
%   wITPC is a little more convoluted, especially since it can't be used
%   with single-trial LME like all otehr analyses in the paper
%% Run weighted inter-trial phase coherence (wITPC) statistics:
%   wITPC weights single-trial phase by model regressors, which are then
%   z-scored to bootstrapped null distribution within SBJ before
%   t-test against zero at the group level
%   Only for one channel
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load single-trial TFR phase data and design matrix (model regressors, SBJ factor)
%       Optional: z-score model regressors within SBJ
%   Compute wITPC across trials, weighting by model regressor
%   Bootstrap null distribution of wITPC by randomizing regressors
%   Z-score wITPC within SBJ to null distribution
%   Independent samples t-test for zwITPC against zero at group level for
%       each regressor, frequency, and time
%   Correct for multiple comparisons (FDR for regressors, times, frequencies)
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   witpc [float array] - wITPC values for [SBJ, regressor, freq, time]
%   zwitpc [float array] - z-scored wITPC values for [SBJ, regressor, freq, time]
%   qvals [float array] - [n_regressors, n_freq, n_time] p values adjusted for multiple comparisons 
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data 
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~an.complex; error('why run this without ITPC an_vars?'); end

stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.an_style,'wITPC'); error('stat_id not using wITPC!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Data and Build Model
tic
cfgs  = []; cfgs.latency = st.stat_lim;
for s = 1:numel(SBJs)
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    % Load behavior
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    full_cond_idx = fn_condition_index(cond_lab, bhv);
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
    
    % Z-score SBJ model regressors
    model = NaN(size(tmp.model));
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            model(:,reg_ix) = ...
                (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
        end
    elseif st.rank_order
        for reg_ix = 1:numel(reg_lab)
            [~, ~, model(:,reg_ix)] = unique(tmp.model(:,reg_ix));
        end
    else
        model = tmp.model;
    end
    
    % Load phase data
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_' proc_id '_' an_id '.mat'],'tfr');
    if numel(tfr.label)>1; error('assuming single channel for now!'); end
    
    % Select time and trials of interest
    model = model(full_cond_idx,:);
    cfgs.trials  = find(full_cond_idx);
    st_tfr = ft_selectdata(cfgs, tfr);
    
    % Initialize matrices now that we know time-frequency axes
    if s==1
        time_vec   = st_tfr.time;
        fois       = st_tfr.freq;
        witpc      = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec)]);
        zwitpc     = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec)]);
        witpc_null = nan([numel(SBJs) numel(reg_lab) numel(fois) numel(time_vec) st.n_boots]);
    end
    
    % Weighted Inter-Trial Phase Coherence
    trial_idx_perm = zeros([size(model,1) st.n_boots]);
    for reg_ix = 1:numel(reg_lab)
        fprintf('%s (%d/%d) freq: ',reg_lab{reg_ix},reg_ix,numel(reg_lab));
        for f_ix = 1:numel(fois)
            fprintf('%.3f..',fois(f_ix));
            for t_ix = 1:numel(time_vec)
                % Compute weighted ITPC across trials (within SBJ)
                %   ITPC is the magnitude of the mean phase angle vector after
                %   using Euler's formula to bring phase angles into polar coordinate space
                % Before averaging into ITPC, weight single-trial phase vectors
                %   by model regressor
                witpc(s,reg_ix,f_ix,t_ix) = ...
                    abs(mean(exp(1i*angle(st_tfr.fourierspctrm(:,1,f_ix,t_ix))).*...
                    model(:,reg_ix)));
                
                % Generate null distribution
                %   wITPC values depend on the model regressor, so they
                %   should be normalized (z-score) before group statistics
                for boot_ix = 1:st.n_boots
                    % Only define randomization once per SBJ
                    if reg_ix==1 && f_ix==1 && t_ix==1
                        trial_idx_perm(:,boot_ix) = randperm(size(model,1));
                    end
                    % Re-compute wITPC with randomized regressors
                    witpc_null(s,reg_ix,f_ix,t_ix,boot_ix) = ...
                        abs(mean(exp(1i*angle(st_tfr.fourierspctrm(:,1,f_ix,t_ix))).*...
                        model(trial_idx_perm(:,boot_ix),reg_ix)));
                end
                % Z-score wITPC values
                zwitpc(s,reg_ix,f_ix,t_ix) = ...
                    (witpc(s,reg_ix,f_ix,t_ix)-mean(witpc_null(s,reg_ix,f_ix,t_ix,:),5)) / ...
                    std(witpc_null(s,reg_ix,f_ix,t_ix,:),[],5);
            end
        end
        fprintf('\n');
    end
    fprintf('End of SBJ %d / %d: %s\n',s,numel(SBJs),SBJs{s});
    toc
    
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
            % Compare group z-scored wITPC to zero
            [~, pvals(reg_ix,f_ix,t_ix)] = ttest(squeeze(zwitpc(:,reg_ix,f_ix,t_ix)));
        end
    end
end

% Correct for Multiple Comparisons (regressors, times, frequencies)
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
stat_out_fname = [stat_out_dir SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','witpc','zwitpc','qvals','SBJs');
% zwitpc_null is too big and not worth saving, just the zwitpc

end
