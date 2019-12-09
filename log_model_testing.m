SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
           'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};
proc_id = 'eeg_full_ft';

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%%
addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data 
wl_fig = figure;
tol_fig = figure;
[n_rc,~] = fn_num_subplots(numel(SBJs));
betas = nan([numel(SBJs) 2]);
total_score = nan(size(SBJs));
for s = 1:numel(SBJs)
    SBJ = SBJs{s};
%     eval(['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m']);
    
    % Load Data
    load([root_dir 'PRJ_Error_eeg/data/' SBJ '/03_events/' SBJ '_behav_' proc_id '_final.mat']);
    total_score(s) = sum(bhv.score);
    
    % Select Data
    s_idx = fn_condition_index({'Su'},bhv);
    X = bhv.tol(~s_idx);
    y = double(bhv.hit(~s_idx));
    [~,sort_idx] = sort(X);
    
    % Logistic regression
    betas(s,:) = glmfit(X,y,'binomial','link','logit');
    
    z = betas(s,1) + (X * betas(s,2));
    log_z = 1 ./ (1+exp(-z));
    
    % Plot Results
    set(0, 'currentfigure', wl_fig);
    subplot(n_rc(1),n_rc(2),s); hold on;
    scatter(X, y);
    plot(X(sort_idx),log_z(sort_idx));
    xlim([0 0.4]);
    ylim([-0.1 1.1]);
    xlabel('Tolerance (s)');
    ylabel('Loss/Hit');
    title([SBJ ': beta=' num2str(betas(s,2),'%.1f')]);
    
    % Plot model vs. tolerance
    set(0, 'currentfigure', tol_fig);
    subplot(n_rc(1),n_rc(2),s); hold on;
    scatter(X,log_z);
    plot(X(sort_idx),log_z(sort_idx));
    xlabel('Tolerance (s)');
    ylabel('p(Win)');
    title([SBJ ': beta=' num2str(betas(s,2),'%.1f')]);
end

[tmp,p] = corrcoef(total_score,betas(:,2));
figure;
scatter(total_score,betas(:,2));
xlabel('Total Score');
ylabel('beta');
title(['r = ' num2str(tmp(1,2),'%.3f') ' (p = ' num2str(p(1,2),'%.3f') ')']);

