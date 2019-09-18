function CALC_Individ_FRNpeaks(SBJ, proc_id, plt_id, an_id)
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
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stats_fname = [SBJ_vars.dirs.proc SBJ '_' an_id 'stats_4cond.mat'];
load(stats_fname);
fig_name = [SBJ '_4Conds_Cz'];
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' SBJ '/4_Conditions/' an_id '/'];
for x = 1:numel(roi_erp)
    openfig([fig_dir fig_name '.fig'])
    findpeaks(roi_erp{x}.avg) %need this to get arrows
    [peak, index] = findpeaks(roi_erp{x}.avg);
    disp(numel(peak));
    index_pos(x) = input('Please type in the index of the peak of the first positivity!');
    close
    peak_pos(x) = peak(index_pos(x));
    openfig([fig_dir fig_name '.fig'])
    findpeaks(-roi_erp{x}.avg) %need this to get arrows
    [peak, index] = findpeaks(-roi_erp{x}.avg);
    disp(numel(peak));
    index_neg(x) =  input('Please type in the index of the peak of the following negativity!');
    close
    peak_neg(x) = -peak(index_neg(x));
    peak_diff(x) = peak_pos(x) - peak_neg(x);
end   
stats_fname = [SBJ_vars.dirs.proc SBJ '_' an_id 'FRNdiffs_4cond.mat'];
save(stats_fname, '-v7.3', 'peak_diff');
 