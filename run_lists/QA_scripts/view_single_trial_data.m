%% View Single Trial Contributions to ERPs
SBJ = 'EEG04';
erp_id = 'ERP_Fz_F2t1_dm2t0_fl05t20';
pow_id = 'POW_Fz_F2t1_dm2t0_fl4t8';
proc_id = 'eeg_full_ft';

load([root_dir 'PRJ_Error_eeg/data/' SBJ '/04_proc/' SBJ '_' an_id '.mat']);
load([root_dir 'PRJ_Error_eeg/data/' SBJ '/03_events/'  SBJ '_behav_' proc_id '_final.mat']);

%% Select conditions (and trials)
[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(cond_lab, bhv);

cond_mat = [cond_idx; 1-bhv.rt];
trials = nan([numel(roi.trial) length(roi.time{1})]);
for t = 1:numel(roi.trial)
    trials(t,:) = roi.trial{t};
end

%%

plot(roi.time{1},trials);hold on; plot(roi.time{1},mean(trials),'color','k','linewidth',3)

figure;imagesc(roi.time{1},1:size(trials,1),trials)
set(gca,'YDir','normal')