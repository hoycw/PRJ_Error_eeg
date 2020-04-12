function group_accuracy_plots_TT(SBJ_id, proc_id)
% Plots group accuracy for TT Task as a bar plot
%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Necessary Data
sbj_list = load_SBJ_file(SBJ_id);
SBJ_colors = distinguishable_colors(numel(sbj_list));
cond_lab = {'easy','hard'};
acc_cond = zeros([numel(sbj_list) 2]);

for sbj_ix = 1:numel(sbj_list)
    SBJ = sbj_list{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load Behavior Data
    if any(strcmp(SBJ, {'EP01','EP02','EP03','EP04','EP05'}))
        [bhv] = fn_load_behav_csv_old([SBJ_vars.dirs.events SBJ '_behav.csv']);
    else
        [bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);
    end
    % Calculate Subject Level Accuracy
    for cond_ix = 1:2
        acc_cond(sbj_ix,cond_ix) = mean(bhv.hit(strcmp(bhv.cond,cond_lab{cond_ix})));
    end
end

%% Calculate Average across Subjects
for cond_ix = 1:2
    avg_cond(cond_ix) = mean(acc_cond(:,cond_ix));
    std_dev_cond(cond_ix) = std(acc_cond(:,cond_ix));
end
%% Plot Accuracy per Condition Across Subjects
figure;
labels = categorical({'Easy', 'Hard'});
bar(labels,avg_cond)                
hold on
er = errorbar(1:2,avg_cond,std_dev_cond,std_dev_cond);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
fig_name = [root_dir 'PRJ_Error_eeg/results/BHV/accuracy/GRP_accuracy_TT.png'];
fprintf('Saving %s\n',fig_name);
saveas(gcf,fig_name);

%% Plot Accuracy per Subject
figure; hold on;
for sbj_ix = 1:numel(sbj_list)
    line([1 2],acc_cond(sbj_ix,:),'Color',SBJ_colors(sbj_ix,:));
    scatter([1 2], acc_cond(sbj_ix,:),30,SBJ_colors(sbj_ix,:));
end
ax = gca;
ax.XLabel.String = 'Condition';
ax.XTick = [1 2];
ax.XTickLabels = cond_lab;
ax.XLim = [0.5 2.5];

ax.YLim = [0 1];
ax.YLabel.String = 'Accuracy';

set(ax,'FontSize',16');

fig_name = [root_dir 'PRJ_Error_eeg/results/BHV/accuracy/GRP_plot_SBJ_accuracy_TT.png'];
fprintf('Saving %s\n',fig_name);
saveas(gcf,fig_name);
