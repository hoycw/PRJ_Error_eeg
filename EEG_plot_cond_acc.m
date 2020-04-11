eeg_behav_SBJ = {'EEG01','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12',...%'EEG02',
    'EP07','EP09','EP14','EP15','EP16','EP17','EP18','EP19'};%'EP08','EP10','EP11',
proc_id = 'eeg_full_ft';
sbj_list = eeg_behav_SBJ;
SBJ_colors = distinguishable_colors(numel(sbj_list));
cond_lab = {'easy','hard'};
acc_cond = zeros([numel(sbj_list) 2]);
for sbj_ix = 1:numel(sbj_list)
SBJ = sbj_list{sbj_ix};
fprintf('================= Processing: %s =================\n',SBJ);
% Load variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

% Compute mean RT
%     data_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
%     load(data_fname,'ignore_trials');
if any(strcmp(SBJ, {'EP01','EP02','EP03','EP04','EP05'}))
[bhv] = fn_load_behav_csv_old([SBJ_vars.dirs.events SBJ '_behav.csv'], []);
else
[bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv'], []);
end
for cond_ix = 1:2
acc_cond(sbj_ix,cond_ix) = mean(bhv.hit(strcmp(bhv.cond,cond_lab{cond_ix})));
end
end

%%
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
