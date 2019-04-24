function trl = tt_trialfun(cfg);
% Trial cutting function for Target Time EEG experiment
% Required Fields:
%   cfg.dataset [str] - pathname to dataset from which to read the events
%   cfg.trialdef.eventtype [str] - name of the channel with the event codes
%   cfg.trialdef.eventvalue [int] - code for desired segmenting events
%   cfg.trialdef.prestim [float] - latency in sec to cut pre-event
%   cfg.trialdef.poststim [float] - latency in sec to cut post-event
% Optional Fields:
%   cfg.resamp_freq [int] - sampling rate if data was resampled (e.g., downsampling)

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset, 'header', hdr);
fprintf('%i events found!\n',numel(event));

% Compute resampling factor
if isfield(cfg,'resamp_freq')
    resamp_factor = cfg.resamp_freq/hdr.Fs;
    srate = cfg.resamp_freq;
else
    resamp_factor = 1;
    srate = hdr.Fs;
end

% Find event cuts and built trl matrix
trl = [];
for i=1:length(event)
  if strcmp(event(i).type, cfg.trialdef.eventtype)
    % it is a trigger, see whether it has the right value
    if ismember(event(i).value, cfg.trialdef.eventvalue)
      % add this to the trl definition
      begsample     = event(i).sample*resamp_factor + round(cfg.trialdef.prestim*srate);
      endsample     = event(i).sample*resamp_factor + round(cfg.trialdef.poststim*srate)-1;
      offset        = cfg.trialdef.prestim*srate;  
      trigger       = event(i).value; % remember the trigger (=condition) for each trial
      trl(end+1, :) = [round([begsample endsample offset])  trigger]; 
    end
  end
end
