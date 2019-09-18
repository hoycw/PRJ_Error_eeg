function CALCPLOT_Grp_FRNpeaks(SBJs, proc_id, plt_id, an_id)
%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);
addpath(ft_dir);
ft_defaults
%% Load preprocessed data
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
for s = 1:length(SBJs)
    SBJ_vars_cmd{s} = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd{s});
    SBJ_vars_all{s} = SBJ_vars;
end
for s = 1:length(SBJs)
    clean_stats_fname{s} = [SBJ_vars_all{s}.dirs.proc SBJs{s} '_' an_id 'FRNdiffs_4cond.mat'];
end
peak_diffs= cell(size(SBJs));
for s = 1:length(SBJs)
    tmp = load(clean_stats_fname{s}); peak_diffs{s} = tmp.peak_diff;
end
for x = 1:4
    peak_diff_total = 0;
    for s = 1:numel(SBJs)
        peak_diff_conds(s,x) = peak_diffs{s}(x);
        peak_diff_total = peak_diffs{s}(x)+peak_diff_total;
    end
    peak_diff_mean(x) = peak_diff_total/numel(SBJs);
end
x=categorical({'Correct/Easy', 'Error/Easy', 'Correct/Hard', 'Error/Hard'});
x = reordercats(x,{'Correct/Easy', 'Error/Easy', 'Correct/Hard', 'Error/Hard'});
bar(x, peak_diff_mean)
set(gca,'xticklabel',{'Correct/Easy', 'Error/Easy', 'Correct/Hard', 'Error/Hard'})
ylabel('Average Microvolts');

[h(1),p(1),ci,stats] = ttest(peak_diff_conds(:,1), peak_diff_conds(:, 2));
[h(2),p(2),ci,stats] = ttest(peak_diff_conds(:,1), peak_diff_conds(:, 3));
[h(3),p(3),ci,stats] = ttest(peak_diff_conds(:,1), peak_diff_conds(:, 4));
[h(4),p(4),ci,stats] = ttest(peak_diff_conds(:,2), peak_diff_conds(:, 3));
[h(5),p(5),ci,stats] = ttest(peak_diff_conds(:,2), peak_diff_conds(:, 4));
[h(6),p(6),ci,stats] = ttest(peak_diff_conds(:,3), peak_diff_conds(:, 4));

index = find(h == 1);
if index
    fprintf(index);
end
end
    
    
    
    %%
%{
for cond_ix = 1:numel(conds)
    for sbj_ix = 1:length(SBJs)
        peak_diff
    end
%}
