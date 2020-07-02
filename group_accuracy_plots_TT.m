function group_accuracy_plots_TT(SBJ_id, conditions, fig_ftype)
%% Plots group accuracy for Target Time task as a bar plot with SBJ lines overlay

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
SBJs = fn_load_SBJ_list(SBJ_id);
SBJ_colors = distinguishable_colors(numel(SBJs));
[cond_lab, cond_names, ~, ~, ~] = fn_condition_label_styles(conditions);
acc_cond = zeros([numel(SBJs) 2]);

for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    fprintf('================= Processing: %s =================\n',SBJ);
    
    % Load Behavior Data
    if any(strcmp(SBJ, {'EP01','EP02','EP03','EP04','EP05'}))
        error('Why run with EP01,02,03,04,05 (original pilot SBJ)?');
        % [bhv] = fn_load_behav_csv_old([root_dir 'PRJ_Error_eeg/data/' SBJ '/03_events/' SBJ '_behav.csv']);
    else
        [bhv] = fn_load_behav_csv([root_dir 'PRJ_Error_eeg/data/' SBJ '/03_events/' SBJ '_behav.csv']);
    end
    
    % Calculate Subject Level Accuracy
    s_idx = fn_condition_index({'Su'},bhv);
    training_idx = bhv.blk==0;
    
    for cond_ix = 1:2
        cond_idx = all([fn_condition_index(cond_lab(cond_ix),bhv) ~s_idx ~training_idx],2);
        acc_cond(sbj_ix,cond_ix) = mean(bhv.hit(cond_idx));
    end
end

%% Calculate Average across Subjects
avg_cond = nan(size(cond_lab));
std_cond = nan(size(cond_lab));
% sem_cond = nan(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    avg_cond(cond_ix) = mean(acc_cond(:,cond_ix));
    std_cond(cond_ix) = std(acc_cond(:,cond_ix));
    % sem_cond(cond_ix) = std_cond(cond_ix)/sqrt(numel(SBJs));
    fprintf('%s %s mean +/- SD = %.04f +/- %.04f\n',SBJ_id,cond_lab{cond_ix},...
        avg_cond(cond_ix),std_cond(cond_ix));
end

%% Plot Accuracy per Condition Across Subjects
fig_dir = [root_dir 'PRJ_Error_eeg/results/BHV/accuracy/'];
fig_name = [SBJ_id '_accuracy_' conditions];
figure('Name',fig_name); hold on

% Plot Bars
bar(1:numel(cond_lab),avg_cond);

% Plot single SBJ data
for sbj_ix = 1:numel(SBJs)
    line(1:numel(cond_lab), acc_cond(sbj_ix,:), 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
end

% Plot error bars
er = errorbar(1:numel(cond_lab),avg_cond,std_cond./2);
er.Color = [0 0 0];
er.LineStyle = 'none';

set(gca,'XLim',[0.5 numel(cond_lab)+0.5]);
set(gca,'XTick',1:numel(cond_lab));
set(gca,'XTickLabel',cond_names);
set(gca,'YLim',[0 1]);
ylabel('Accuracy');
title([SBJ_id ' Accuracy Across Conditions']);
set(gca,'FontSize',16);

fig_fname = [fig_dir fig_name '.' fig_ftype];
fprintf('Saving %s\n',fig_fname);
saveas(gcf,fig_fname);

%% Plot Accuracy per Subject
% figure; hold on;
% for sbj_ix = 1:numel(SBJs)
%     line([1 2],acc_cond(sbj_ix,:),'Color',SBJ_colors(sbj_ix,:));
%     scatter([1 2], acc_cond(sbj_ix,:),30,SBJ_colors(sbj_ix,:));
% end
% ax = gca;
% ax.XLabel.String = 'Condition';
% ax.XTick = [1 2];
% ax.XTickLabels = cond_lab;
% ax.XLim = [0.5 2.5];
% 
% ax.YLim = [0 1];
% ax.YLabel.String = 'Accuracy';
% 
% set(ax,'FontSize',16');
% 
% fig_name = [root_dir 'PRJ_Error_eeg/results/BHV/accuracy/GRP_plot_SBJ_accuracy_TT.' fig_ftype];
% fprintf('Saving %s\n',fig_name);
% %saveas(gcf,fig_name);

end