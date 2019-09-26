function CALC_FRN_Individ(SBJ, proc_id, plt_id, an_id)
%Purpose: This function loads the data for each condition for a given subject.  Then it allows users to input the index of the peak that represents the index of the first positivity and then the corresponding negativity.  Then it computes the difference between the two values and stores them.
%Inputs
%SBJ = string (ie) 'EEG01'
%proc_id = 'eeg_full_ft'
%plt_id = 'ts_F15to28_evnts_sigPatch'
%an_id = 'ERP_Cz_F_trl15t28_flt05t20_stat06'
%fig_vis = whether or not the plot should pop up ('on' or 'off')
%save_fig = 0 or 1
%fig_ftype = 'png' (a string)
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
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/4_Conditions/' an_id '/'];
stats_dir = [root_dir 'PRJ_Error_eeg/results/Stats/4_Conditions/' an_id '/'];
stats_fname = [stats_dir an_id SBJ '_4Conds.mat'];
load(stats_fname);
fig_name_FRN = [SBJ '_FRN_' an.ROI{1}];
for x = 1:numel(roi_erp)
    fig_name = [SBJ '_4Conds_' an.ROI{1}];
    fig_path = [fig_dir fig_name '.fig'];
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
stats_fname = [stats_dir an_id SBJ '_FRN.mat'];
save(stats_fname, '-v7.3', 'peak_diff');
 
