function [trials] = fn_ft_cut_trials_equal_len(data,events,event_types,lim)
%% Cut continuous data into trials using ft_redefinetrial
%   This version only takes one set of events and cuts all trials to equal
%   lengths (e.g., stimulus or response locked, but taking only one into account).
%
% Inputs:
%   data [ft data struct]- continuous data in fieldtrip format
%   events [int array]- array of sample indices to lock cutting
%   event_types [int array]- array of condition labels as integers
%   lim [int tuple array]- signed, two integer array specifying how to cut around events
%           lim(1) = baseline, e.g. -500 would be 500 samples before events
%           lim(2) = post-event length, e.g. 1000 would be 1000 samples after events
% Outputs:
%   trials [ft data struct]- a segmented fieldtrip data structure
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath(ft_dir);
ft_defaults

% Check size of inputs
if numel(events)~=numel(event_types)
    error('ERROR: length of events and event_types mismatch');
end
if numel(lim)~=2
    error('ERROR: lim must have only 2 elements- pre- and post-event lengths');
end

% Cut into trials
cfg = [];
% cfg.dataset = strcat(SBJ_dir,'04_proc/',data_id,'_data_only.mat');
% Define trials as my whole period of interest (buffer to buffer, no offset)
%   These numbers are in samples!
cfg.trl = [events+lim(1),...        % start of the trial including baseline
           events+lim(2),...    % end of trial
           repmat(lim(1),length(events),1),...         %time of event relative to start of trial
           event_types];                 % trial type
% cfg.continuous = 'yes';
trials = ft_redefinetrial(cfg, data);


end
