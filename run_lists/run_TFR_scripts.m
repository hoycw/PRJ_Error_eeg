function run_TFR_scripts
%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%%
addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% General parameters
SBJs = {'EP06','EP07','EP08','EP09','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};

%% Run preprocessing
proc_id = 'eeg_full_ft';
an_id = 'TFR_Cz_F2t1_z2t0_fl2t40';
save_fig    = 1;
fig_vis     = 'off';


for s = 1:numel(SBJs)
     SBJ05a_TFR_save(SBJs{s}, proc_id, an_id)
     SBJ05b_TFR_plot(SBJs{s}, 'DifOut', proc_id, an_id, save_fig);
     SBJ05b_TFR_plot(SBJs{s}, 'Out', proc_id, an_id, save_fig);
end
end