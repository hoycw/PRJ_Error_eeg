function SBJ04b_RL_model_plot_grp_fits(SBJs,proc_id,stat_id)
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
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, cond_colors, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

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

%% Load Data and Build Model
model = nan([sum(n_trials) numel(reg_lab)]);
for r_ix = 1:numel(reg_lab)
    model(:,r_ix) = eval(reg_lab{r_ix});
end

%% Compute and plot correlations between regressors
reg_corr = corr(model,'rows','complete');

% Plot design matrix
fig_name = [SBJ '_' model_id '_design'];
figure('Name',fig_name);
imagesc(model);
xticklabels(reg_lab);
colorbar;
saveas(gcf,[SBJ_vars.dirs.proc fig_name '.png']);

% Plot regressor correlation matrix
fig_name = [SBJ '_' model_id '_design_corr'];
figure('Name',fig_name);
imagesc(reg_corr);
xticklabels(reg_lab);
yticklabels(reg_lab);
colorbar;
saveas(gcf,[SBJ_vars.dirs.proc fig_name '.png']);

%% Plot Regressors by Condition
cond_idx = fn_condition_index(cond_lab, bhv);

fig_name = [SBJ '_' model_id '_reg_cond'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 1 1]);
[n_rc,~] = fn_num_subplots(numel(reg_lab));

for reg_ix = 1:numel(reg_lab)
    subplot(n_rc(1),n_rc(2),reg_ix); hold on;
    violin_data = struct;
    for cond_ix = 1:numel(cond_lab)
        violin_data.(cond_lab{cond_ix}) = model(cond_idx==cond_ix,reg_ix);
    end
    violins = violinplot(violin_data);
    
    % Plot properties
    legend_obj = cell(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        % Change scatter colors to mark condition
        violins(cond_ix).ViolinColor = cond_colors{cond_ix};
        violins(cond_ix).BoxPlot.FaceColor = cond_colors{cond_ix};
        violins(cond_ix).EdgeColor = cond_colors{cond_ix};
        % Grab violin for legend
        legend_obj{cond_ix} = violins(cond_ix).ViolinPlot;
    end
    title(reg_lab{reg_ix});
    if strcmp(reg_lab{reg_ix},'pWin'); legend([legend_obj{:}],cond_lab); end
    set(gca,'FontSize',16);
end
saveas(gcf,[SBJ_vars.dirs.proc fig_name '.png']);

%% Save Results
stat_out_fname = [SBJ_vars.dirs.proc SBJ '_model_' model_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','model');

end
