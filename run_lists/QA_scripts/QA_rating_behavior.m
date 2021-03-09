%% Target Time Performance Rankings
root_dir='/Volumes/hoycw_clust/';
SBJ_id = 'ratings_all';
proc_id = 'eeg_full_ft';

SBJs = fn_load_SBJ_list(SBJ_id);

%% Convert labels
SBJ_lab = SBJs;


%%
addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load data
tol = nan([numel(SBJs) 3]);
acc = nan([numel(SBJs) 3]);
pts = nan([numel(SBJs) 3]);
out = nan([numel(SBJs) 1]);
for s = 1:numel(SBJs)
    % Load data
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    bhvs{s} = tmp.bhv;
    
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_' proc_id '_rate01_orig_exclude_trial_ix.mat'],'rt_low_ix','rt_high_ix');
    out(s) = sum([numel(rt_low_ix) numel(rt_high_ix)]);
    
    ez_idx = strcmp(bhvs{s}.cond,'easy');
    
    tol(s,1) = nanmean(bhvs{s}.tol(ez_idx));
    tol(s,2) = nanmean(bhvs{s}.tol(~ez_idx));
    tol(s,3) = nanmean(bhvs{s}.tol);
    
    acc(s,1) = nanmean(bhvs{s}.hit(ez_idx));
    acc(s,2) = nanmean(bhvs{s}.hit(~ez_idx));
    acc(s,3) = nanmean(bhvs{s}.hit);
    
    pts(s,1) = sum(bhvs{s}.score(ez_idx));
    pts(s,2) = sum(bhvs{s}.score(~ez_idx));
    pts(s,3) = sum(bhvs{s}.score);
end

%% Rank SBJs
[~,tol_sort_idx] = sort(tol(:,3),'ascend');
[~,acc_sort_idx] = sort(acc(:,3),'descend');
[~,pts_sort_idx] = sort(pts(:,3),'descend');
[~,out_sort_idx] = sort(out,'descend');

%% Print results
fprintf('Best Tolerance:\n');
for s = 1:numel(SBJs)
    fprintf('\t(%i) %s: %.3f (%.3f easy, %.3f hard)\n',s,SBJ_lab{tol_sort_idx(s)},...
        tol(tol_sort_idx(s),3),tol(tol_sort_idx(s),1),tol(tol_sort_idx(s),2));
end

fprintf('\n\n');
fprintf('Best Accuracy:\n');
for s = 1:numel(SBJs)
    fprintf('\t(%i) %s: %.3f (%.3f easy, %.3f hard)\n',s,SBJ_lab{acc_sort_idx(s)},...
        acc(acc_sort_idx(s),3),acc(acc_sort_idx(s),1),acc(acc_sort_idx(s),2));
end

fprintf('\n\n');
fprintf('Best Score:\n');
for s = 1:numel(SBJs)
    fprintf('\t(%i) %s: %.3f (%.3f easy, %.3f hard)\n',s,SBJ_lab{pts_sort_idx(s)},...
        pts(pts_sort_idx(s),3),pts(pts_sort_idx(s),1),pts(pts_sort_idx(s),2));
end

fprintf('\n\n');
fprintf('Worst # Space Outs:\n');
for s = 1:numel(SBJs)
    fprintf('\t(%i) %s: %.3f\n',s,SBJ_lab{out_sort_idx(s)},out(out_sort_idx(s)));
end

%% Plot Outliers
out_thresh3 = mean(out) + std(out)*3;
out_thresh2 = mean(out) + std(out)*2;
figure;
histogram(out,25);
xlabel('# RT outliers');
ylabel('SBJ count');
avg_line = line([mean(out) mean(out)],ylim,'Color','k','LineWidth',2);
thresh2_line = line([out_thresh2 out_thresh2],ylim,'Color','r','LineWidth',2);
thresh3_line = line([out_thresh3 out_thresh3],ylim,'Color','r','LineWidth',2);
legend([avg_line thresh2_line thresh3_line],{['Mean=' num2str(mean(out))],...
    ['2SDs=' num2str(out_thresh2)],['3SDs=' num2str(out_thresh3)]});
set(gca,'FontSize',16);

%% Plot results
figure('units','normalized','OuterPosition',[0 0 1 1]);
subplot(4,1,1); hold on;
bar(tol(tol_sort_idx,3));
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJ_lab(tol_sort_idx));
ylabel('Average Target Zone Size (sec)');
title('Best Overall Skill (Smallest Target)');
set(gca,'FontSize',16);

subplot(4,1,2); hold on;
acc_fudge = 0.01;
bar(acc(acc_sort_idx,3));
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJ_lab(acc_sort_idx));
ylim([min(acc(:,3))-acc_fudge max(acc(:,3))+acc_fudge]);
ylabel('Total Accuracy');
title('Best Total Accuracy');
set(gca,'FontSize',16);

subplot(4,1,3); hold on;
pts_fudge = 500;
bar(pts(pts_sort_idx,3));
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJ_lab(pts_sort_idx));
ylim([min(pts(:,3))-pts_fudge max(pts(:,3))+pts_fudge]);
ylabel('Total Score');
title('Best Total Score');
set(gca,'FontSize',16);

subplot(4,1,4); hold on;
bar(out(out_sort_idx));
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJ_lab(out_sort_idx));
ylabel('# Spaced Out Trials');
title('Most Spaced Out (Reaction Times Off by > 0.4 sec)');
set(gca,'FontSize',16);
