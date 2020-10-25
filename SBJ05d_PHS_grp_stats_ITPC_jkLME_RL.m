function SBJ05d_PHS_grp_stats_ITPC_jkLME_RL(SBJ_id,proc_id,an_id,stat_id)
%% Run Linear Mixed-Effects multiple regression model on single-trial 
% inter-trial phase clustering (ITPC) at each time-frequency point for across all SBJs and trials
%   Uses jackknife procedure to obtain single-trial contributions to ITPC
%   Only for one channel
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load single-trial TFR phase data and design matrix (model regressors, SBJ factor)
%       Optional: z-score model regressors within SBJ
%   Compute ITPC across trials, then jackknife single-trial contributions
%   Run LME RL model per time-frequency point
%   Correct for multiple comparisons (FDR for regressors, times, frequencies)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class for each time-frequency point
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
if ~strcmp(st.an_style,'lme'); error('stat_id not using lme!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(st.stat_cond);

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
tic
cfgs  = []; cfgs.latency = st.stat_lim;
model = zeros([sum(n_trials) numel(reg_lab)]);
sbj_factor  = zeros([sum(n_trials) 1]);
for s = 1:numel(SBJs)
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    % Load TFR
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_' proc_id '_' an_id '.mat'],'tfr');
    if numel(tfr.label)>1; error('assuming single channel for now!'); end
    
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_tfr = ft_selectdata(cfgs, tfr);
    
    % Get SBJ index and initialize matrices now that we know time axis
    if s==1
        time_vec = st_tfr.time;
        fois     = st_tfr.freq;
        if strcmp(st.measure,'ts')
            data  = nan([sum(n_trials) numel(fois) numel(time_vec)]);
        elseif strcmp(st.measure,'mean')
            error('st.measure should not be mean for phase angles!');
        else; error(['unknown st.measure: ' st.measure]);
        end
        
        sbj_idx = 1:n_trials(s);
    else
        sbj_idx = sum(n_trials(1:s-1))+1:sum(n_trials(1:s));
    end
    
    % Track SBJ in design matrix
    sbj_factor(sbj_idx) = s*ones([n_trials(s) 1]);
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' st.model_id '.mat']);
    
    % Z-score SBJ model regressors
    sbj_model = NaN([sum(full_cond_idx{s}~=0) size(tmp.model,2)]);
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            sbj_model(:,reg_ix) = ...
                (tmp.model(full_cond_idx{s}~=0,reg_ix)-nanmean(tmp.model(full_cond_idx{s}~=0,reg_ix)))./...
                nanstd(tmp.model(full_cond_idx{s}~=0,reg_ix));
        end
    else
        sbj_model = tmp.model(full_cond_idx{s}~=0,:);
    end
    model(sbj_idx,:) = sbj_model;
    
    % Jack-Knife Inter-Trial Phase Coherence
    if strcmp(st.measure,'ts')
        fprintf('%s freq: ',SBJs{s});
        
        itpc = nan([numel(fois) numel(time_vec)]);
        for f_ix = 1:numel(fois)
            fprintf('%.3f..',fois(f_ix));
            for t_ix = 1:numel(time_vec)
                % Compute ITPC across trials (within SBJ)
                %   ITPC is the magnitude of the mean phase angle vector after
                %   using Euler's formula to bring phase angles into polar coordinate space
                itpc(f_ix,t_ix) = abs(mean(exp(1i*angle(st_tfr.fourierspctrm(:,1,f_ix,t_ix)))));
                
                % Jack-knife ITPC for each trial contribution
                %   Recompute ITPC after excluding one trial, then subtract
                %   from ITPC with all trials to get single-trial contribution
                for trl_ix = 1:n_trials(s)
                    jk_idx = setdiff(1:n_trials(s),trl_ix);
                    jk_itpc = abs(mean(exp(1i*angle(st_tfr.fourierspctrm(jk_idx,1,f_ix,t_ix)))));
                    data(sbj_idx(trl_ix),f_ix,t_ix) = itpc(f_ix,t_ix) - jk_itpc;
                end
            end
        end
        fprintf('\n');
    elseif strcmp(st.measure,'mean')
        error('st.measure should not be mean for phase angles!');
    else; error(['unknown st.measure: ' st.measure]);
    end
    fprintf('End of %s (%d / %d): %s\n',SBJs{s},s,numel(SBJs));
    toc
    
    clear tmp tfr st_tfr sbj_model itpc jk_idx jk_itpc
end
fprintf('\t\t SBJ jack-knife ITPC complete:');
toc

%% Build table
fprintf('========================== Running Stats ==========================\n');
tic

% Build Model Table
tbl = table;
for reg_ix = 1:numel(reg_lab)
    tbl.(reg_lab{reg_ix}) = model(:,reg_ix);
end
tbl.SBJ = categorical(sbj_factor);

% Create Model Formula
reg_formula = strjoin(reg_lab,' + ');
formula = ['PHS ~ ' reg_formula ' + (1|SBJ)'];  % random intercepts for SBJ

% Run Model
lme = cell([numel(fois) numel(time_vec)]);
pvals = nan([numel(reg_lab) numel(fois) numel(time_vec)]);
for f_ix = 1:numel(fois)
    for t_ix = 1:numel(time_vec)
        tbl.PHS = data(:,f_ix,t_ix);
        lme{f_ix,t_ix} = fitlme(tbl,formula);
        pvals(:,f_ix,t_ix) = lme{f_ix,t_ix}.Coefficients.pValue(2:end);
    end
end

% Correct for Multiple Comparisons (regressors, times, frequencies)
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2)*size(pvals,3) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2) size(pvals,3)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

fprintf('\t\t Stats Complete:');
toc

%% Save Results
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(stat_out_dir,'dir')
    [~] = mkdir(stat_out_dir);
end
stat_out_fname = [stat_out_dir SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals','SBJs');

end
