function overlap_trial_ix = fn_find_trials_overlap_epochs(epochs,sample_idx,events,trial_lim)
%% Find trials that overlap with epochs by comparing sample indices
%   Mainly used to find trials that overlap with visually identified artifacts
% INPUTS:
%   epochs [Nx2 int array] - array of [start,stop] pairs of sample indicies marking the epochs
%   sample_idx [Nx1 int array] - array of sample indicies for the
%       continuous data, to be cut into trials based on events+trial_lim
%       e.g., 1:size(data.trial{1},2)
%   events [int array] - sample indices indicating the segmentation point for trials
%   trial_lim [Nx2 int array] - # data points to include [before, after] the event
% OUTPUTS:
%   overlap_trials [Nx1 int array] - column vector of trial indices that overlap with any epochs
[root_dir,~] = fn_get_root_dir();
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);

% Compile list of samples covered by epochs
epoch_samples = [];
for epoch_ix = 1:size(epochs,1)
    epoch_samples = [epoch_samples epochs(epoch_ix,1):epochs(epoch_ix,2)];
end

% Convert trials to sample indices
trials_sample_idx = fn_epoch_cuts_datapad(sample_idx,events,events,trial_lim);

% Find trials that overlap with any epoch indices
overlap_trial_ix = [];
for trial_ix = 1:size(trials_sample_idx,1)
    if ~isempty(intersect(trials_sample_idx(trial_ix,:),epoch_samples))
        overlap_trial_ix = [overlap_trial_ix; trial_ix];
    end
end

end
