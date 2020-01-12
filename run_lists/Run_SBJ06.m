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
SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
           'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};

%% Run Prototype Selection
proc_id   = 'odd_full_ft';
stat_id   = 'CPA';
plt_id    = 'ts_F2to13_evnts_sigLine';

for s = 1: numel(SBJs)
    %ODD02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 'odd_full_ft', 1, 'on', 'ERP_stack_full_events_odd')
    %SBJ02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 1, 'on')
    %SBJ01_preproc(SBJs{x}, 'eeg_full_ft');
    SBJ06a_CPA2(SBJs{s}, proc_id, stat_id, plt_id);
end
