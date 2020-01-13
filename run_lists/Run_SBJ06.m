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
SBJs = {'EEG01','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};
bad_SBJs = {'EEG02'};

%% Run Prototype Selection
proc_id   = 'odd_full_ft';
stat_id   = 'CPA';
plt_id    = 'ts_F2to13_evnts_sigLine';

for s = 1: numel(SBJs)
    %ODD02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 'odd_full_ft', 1, 'on', 'ERP_stack_full_events_odd')
    %SBJ02a_artifact_rejection(SBJs{x}, 'eeg_full_ft', 1, 'on')
    %SBJ01_preproc(SBJs{x}, 'eeg_full_ft');
    SBJ06a_CPA_prototype_selection(SBJs{s}, proc_id, stat_id, plt_id);
end

%% Test Candidate
eeg_proc_id   = 'eeg_full_ft';
odd_proc_id   = 'odd_full_ft';
conditions   = 'DifFB';
an_id     = 'ERP_Fz_F2t1_dm2t0_fl05t20';%,'ERP_Pz_F2t1_dm2t0_fl05t20'};
stat_id   = 'CPA';
plt_id    = 'ts_F2to13_evnts_sigLine';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 1:numel(SBJs)
    SBJ06b_CPA_plot_candidate_ERP(SBJs{s},conditions,eeg_proc_id,odd_proc_id,an_id,stat_id,...
        plt_id,save_fig,'fig_vis',fig_vis,'fig_ftype',fig_ftype);
end
