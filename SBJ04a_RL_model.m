function SBJ04a_RL_model(SBJ,proc_id,stat_id)
% Run reinforcement learning model on single SBJ behavior
% INPUTS:
%   SBJ [str] - ID of subject to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class for each time point

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

[reg_lab, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load and Select Behavior
% Load data
load([root_dir 'PRJ_Error_eeg/data/' SBJ '/03_events/' ...
    SBJ '_behav_' proc_id '_final.mat'],'bhv');
prdm_vars = load([SBJ_vars.dirs.events SBJ '_prdm_vars.mat']);

% Select Conditions of Interest
full_cond_idx = fn_condition_index(cond_lab, bhv);
bhv_fields = fieldnames(bhv);
orig_n_trials = numel(bhv.trl_n);
for f_ix = 1:numel(bhv_fields)
    if numel(bhv.(bhv_fields{f_ix}))==orig_n_trials
        bhv.(bhv_fields{f_ix}) = bhv.(bhv_fields{f_ix})(full_cond_idx~=0);
    end
end
n_trials = numel(bhv.trl_n);
fprintf('\t%s: Loaded %d trials, kept %d for modeling...\n',SBJ,orig_n_trials,n_trials);

%% Run Logistic Regression for Win Prediction
% Select Data (fit on everything except surprise since no outcome)
s_idx = fn_condition_index({'Su'},bhv);
X = bhv.tol(~s_idx);
y = double(bhv.hit(~s_idx));

% Logistic regression
betas = glmfit(X,y,'binomial','link','logit');

z = betas(1) + (bhv.tol * betas(2));
pWin = 1 ./ (1+exp(-z));
expected_score = pWin*2 - 1;
sPE = double(bhv.score)/100 - expected_score;
uPE = abs(sPE);

%% Reaction Time Regressors
% RT (-target_time to center for early/late)
if any(strcmp(reg_lab,'tRT'))
    tRT = bhv.rt-prdm_vars.target;
end

% Previous trial RTs
if any(strcmp(reg_lab,'ptRT'))
    ptRT = nan(size(bhv.rt));
    for t_ix = 2:numel(bhv.rt)
        ptRT(t_ix) = tRT(t_ix-1);
    end
end
if any(strcmp(reg_lab,'p2tRT'))
    p2tRT = nan(size(bhv.rt));
    for t_ix = 3:numel(bhv.rt)
        p2tRT(t_ix) = ptRT(t_ix-2);
    end
end

%% Distance from Tolerance
% Compute distance from closest tolerance bound
%   Positive = win, Negative = loss
early_idx = (bhv.rt-prdm_vars.target)<0;
dist = nan(size(bhv.rt));
dist(early_idx) = bhv.rt(early_idx)-(prdm_vars.target-bhv.tol(early_idx));
dist(~early_idx) = (prdm_vars.target+bhv.tol(~early_idx))-bhv.rt(~early_idx);

% Transform almost hit/miss into large values
if any(strcmp(reg_lab,'lDist'))
    lDist = nan(size(bhv.rt));
    % Separately for wins and losses to maintain correct sign for beta
    % interpretation (+ for win, - for loss)
    lDist(bhv.hit==1) = -log(dist(bhv.hit==1));
    lDist(bhv.hit==0) = log(-dist(bhv.hit==0));
elseif any(strcmp(reg_lab,'iDist'))
    iDist = dist.^-1;
end

% Add unsigned version
if any(strcmp(reg_lab,'ulDist'))
    ulDist = abs(lDist);
elseif any(strcmp(reg_lab,'uiDist'))
    uiDist = abs(iDist);
end

%% Load Data and Build Model
model = nan([sum(n_trials) numel(reg_lab)]);
for r_ix = 1:numel(reg_lab)
    if st.z_reg
        model(:,r_ix) = eval(['(' reg_lab{r_ix} '-nanmean(' reg_lab{r_ix} ...
            '))./nanstd(' reg_lab{r_ix} ');']);
    else
        model(:,r_ix) = eval(reg_lab{r_ix});
    end
end

%% Compute and plot correlations between regressors
reg_corr = corr(model,'rows','complete');

% Plot design matrix
fig_name = [SBJ '_' st.model_lab '_design'];
figure('Name',fig_name);
imagesc(model);
xticklabels(reg_lab);
colorbar;
saveas(gcf,[SBJ_vars.dirs.proc fig_name '.png']);

% Plot regressor correlation matrix
fig_name = [SBJ '_' st.model_lab '_design_corr'];
figure('Name',fig_name);
imagesc(reg_corr);
xticklabels(reg_lab);
yticklabels(reg_lab);
colorbar;
saveas(gcf,[SBJ_vars.dirs.proc fig_name '.png']);

%% Save Results
stat_out_fname = [SBJ_vars.dirs.proc SBJ '_model_' stat_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','model');

end
