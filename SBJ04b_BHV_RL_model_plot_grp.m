function SBJ04b_BHV_RL_model_plot_grp(SBJ_id,proc_id,stat_id, varargin)
% Plots Group Level RL Model Fits (Rereferenced to the midpoint y value and not), as well as the correlation between the Model Fit and several behavioral parameters
%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'save_fig') && ischar(varargin{v+1})
            save_fig = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ~exist('save_fig','var'); save_fig = 1; end

%% Load Data
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

model_id = [st.model_lab '_' st.trial_cond{1}];
[cond_lab, ~, cond_colors, ~, ~] = fn_condition_label_styles(st.trial_cond{1});
%% Load and Select Behavior
% Load data
SBJs = load_SBJ_file(SBJ_id);

% Initialize Variables
sig_x = zeros(numel(SBJs), 401);
sig_y = zeros(numel(SBJs), 401);
new_0_y_ix = zeros(numel(SBJs));
new_0_y_val = zeros(numel(SBJs));
sig_y_aligned = zeros(numel(SBJs), 401);

%% 
for s = 1: numel(SBJs)
    % Load Subject Specific Data
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.proc SBJ '_model_' model_id '.mat']);
    load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'],'bhv');

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
    %% Behavioral Measures
    
    behav_sum.score(s) = sum(bhv.score); % Total Score
    behav_sum.accuracy(s) = sum(bhv.hit)/numel(bhv.hit); % Total Accuracy
    easy_indices = find(strcmp(bhv.cond,'easy')); 
    hard_indices = find(strcmp(bhv.cond,'easy'));
    behav_sum.accuracy_easy(s) = sum(bhv.hit(easy_indices))/numel(easy_indices); % Easy Accuracy
    behav_sum.accuracy_hard(s) = sum(bhv.hit(hard_indices))/numel(hard_indices); % Hard Accuracy

    %% Plot Tolerance vs. Outcome
    sig_step = 0.001;

    % Calculate Model Fit
    sig_x(s,:) = [0:sig_step:0.4];
    sig_y(s,:) = betas(1) + (sig_x(s,:) * betas(2));
    sig_y(s,:) = 1 ./ (1+exp(-sig_y(s,:)));
    
    % Calculate Value for Correlation with Accuracy Analysis -- Maybe also
    % Use beta values?
    y_val_bhv_analysis(s) = max(diff(sig_y(s,:)));
    
    % Calculate Median Centered Values -- find value at midpoint of x (0.2) and
    % align all y values to center around that median value
    %{
    new_0_y_ix(s) = find(sig_x(s,:) == 0.2);
    new_0_y_val(s) = sig_y(s, new_0_y_ix(s));
    sig_y_aligned(s,:) = sig_y(s,:) - new_0_y_val(s);
    %}
    
    % Calculate Median Centered Values -- find value at midpoint of y (0.5) and
    % align all x values to center around that median value
    [~,new_0_x_ix(s)] = min(abs(sig_y(s,:) - 0.5)); % Find x ix for when y = 0.5
    new_0_x_val(s) = sig_x(s, new_0_x_ix(s)); % Find the valeu at that point 
    sig_x_aligned(s,:) = sig_x(s,:) - new_0_x_val(s); % Normalize x values to it

end
%% Plot figs
SBJ_colors = distinguishable_colors(numel(SBJs));
%% Plot All SBJ Level Sigmoids
fig_name = ['GRP_BHV_acc_' model_id '_pWin'];
fig = figure('Name',fig_name,'Visible',fig_vis);
title('GRP Level Accuracy vs. Tolerance');
hold on;
for s = 1:numel(SBJs)
    fit_line = line(sig_x(s,:),sig_y(s,:),'Color', SBJ_colors(s,:));
end
% Figure Parameters
xlabel('Tolerance (s)');
ylabel('Accuracy');
set(gca,'YLim',[0 1]);
set(gca,'FontSize',14);
hold on

%% Normalized Y Values
%{
fig_name_aligned = ['GRP_BHV_acc_' model_id '_pWin_Aligned'];
fig_aligned = figure('Name',fig_name_aligned,'Visible',fig_vis);
title('GRP Level Accuracy vs. Tolerance -- Aligned to Center Point');
hold on;
for s = 1:numel(SBJs)
    fit_line = line(sig_x(s,:),sig_y_aligned(s,:),'Color', SBJ_colors(s,:));
end
% Figure Parameters
xlabel('Tolerance (s)');
ylabel('Accuracy');
set(gca,'FontSize',14);
hold on
%}
%% Normalized X Values
fig_name_aligned = ['GRP_BHV_acc_' model_id '_pWin_Aligned'];
fig_aligned = figure('Name',fig_name_aligned,'Visible',fig_vis);
title('GRP Level Accuracy vs. Tolerance -- Rereferenced');
hold on;
for s = 1:numel(SBJs)
    fit_line = line(sig_x_aligned(s,:),sig_y(s,:),'Color', SBJ_colors(s,:));
end
% Figure Parameters
xlabel('Rereferenced Tolerance (s)');
ylabel('Accuracy');
set(gca,'FontSize',14);
hold on

%% Plot Sigmoid dy/dx vs Behavioral Measures
titles = {'Curve Steepness vs Total Score', 'Curve Steepness vs Total Accuracy', 'Curve Steepness vs Accuracy in Easy Condition'...
    'Curve Steepness vs Accuracy in Hard Condition', };
y_labels = {'Total Score', 'Total Accuracy', 'Accuracy in Easy Condiiton', 'Accuracy in Hard Condition'}
fields = fieldnames(behav_sum);
for ix = 1:numel(fields)
    var_name = fields{ix};
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/model_fits/' model_id '/'];
    fig_name_behav = [fig_dir 'GRP_BHV_acc_' model_id '_pWin_Aligned_' var_name '.' fig_ftype];
    fig_behav = figure('Name',fig_name_behav,'Visible',fig_vis);
    y = getfield(behav_sum, fields{ix});
    line_behav = scatter(y_val_bhv_analysis, y);
    title(titles{ix});
    xlabel('Max Slope of Model Fit');
    ylabel(y_labels{ix});
    fprintf('Saving %s\n',fig_name_behav);
    saveas(fig_behav,fig_name_behav);
end
hold on;
%% Save Figures
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/model_fits/' model_id '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fig_fname_aligned = [fig_dir fig_name_aligned '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname, fig_fname_aligned);
    % Ensure vector graphics if saving
    if any(strcmp(fig_ftype,{'svg','eps'}))
        set(gcf, 'Renderer', 'painters');
    end
    saveas(fig,fig_fname);
    saveas(fig_aligned,fig_fname_aligned);
end
