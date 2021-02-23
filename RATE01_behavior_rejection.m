function RATE01_behavior_rejection(SBJ, proc_id)
%% Reject bad behavior
% (1) Load behavior
% (2) Cut out training, and bad RT trials
% INPUTS:
%   SBJ [str] - name of the SBJ
%   proc_id [str] - name of the preprocessing pipeline parameters (e.g., 'egg_full_ft')
% OUTPUTS: BEHAVIOR
%   bhv [struct] - behavioral data excluding bad trials
%   bhv_fields [cell array] - list of fields in bhv struct
% OUTPUTS: TRIAL EXCLUSION TRACKING
%   NOTE: saving original trial indices means the correct trials can be tossed when
%       loading raw data, e.g., in SBJ05 for TFR filtering on whole data
%   training_ix [array] - original trial indices tossed for training
%   rt_low_ix [array] - original trial indices tossed for too short RTs
%   rt_high_ix [array] - original trial indices tossed for too long RTs
%   exclude_trials [array] - original trial indices tossed for any reason (all unique trial indices combined)

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults
% NOTE: Colin deleted below line (and Sheila created directory/functions) on 8/14/20,
%   hopefully because it's not necessary anymore (didn't run to check)
% addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);

%% Processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

%% Load Behavior
[bhv] = fn_load_rating_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);

%% Exclude bad_trials
% Identify training and bad behavioral trials
training_ix = find(bhv.blk==0);                 % Training block index is 0 here, -1 in python
rt_low_ix   = find(bhv.rt <= proc.rt_bounds(1));
rt_high_ix  = find(bhv.rt >= proc.rt_bounds(2));
exclude_trials = unique(vertcat(training_ix, rt_low_ix, rt_high_ix));
fprintf(2,'\tWarning: Removing %i trials (%i training, %i slow rts, %i fast rts)\n', numel(exclude_trials),...
    numel(training_ix), numel(rt_low_ix), numel(rt_high_ix));

% Exclude bad trials
bhv_fields = fieldnames(bhv);
for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(exclude_trials) = [];
end

%% Save Data
clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv', 'bhv_fields');

% Save excluded trial indices
excluded_fname = [SBJ_vars.dirs.events SBJ '_' proc_id '_rate01_orig_exclude_trial_ix.mat'];
save(excluded_fname,'-v7.3','training_ix','rt_low_ix','rt_high_ix','exclude_trials');

end
