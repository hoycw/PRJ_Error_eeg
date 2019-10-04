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
%checks to see if starting with oddball or not, oddball_section = 1 means yes
oddball_section = 0;
switched = 0;
for i = 1:400
    if (event(i).value == 3)
        oddball_section = 1;
    end
end
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
    if (event(i).value == 255 || event(i).value == 254) && i>4
        % 255 and 254 seem to be interchangable.  Therefore, this checks
        % that it is not the first instance (marking the start of the
        % target_time if you do get a start thing).  HOWEVER, this breaks
        % if you have less than 400 oddballs so this is not a great
        % implementation.  Also this doesn't work for if you start with target_time and it restarts after 400. SS
        % This marks the end of the oddball section and the start of the TT
        oddball_section = abs(oddball_section - 1); %I think this is a cheat way to flip between 0 and 1
        if ~switched
            oddball_section = abs(oddball_section - 1); % I think this is a cheat way to swithc between 0 and 1              
            switched = 1;
        end
    end
    if event(i).value == 254 || event(i).value == 255 && ~oddball_section
         trl = [];
    end
    % Add TT events: event code (1, 2) = (stim, feedback)
    if ~oddball_section && ismember(event(i).value, cfg.trialdef.eventvalue)
      % add this to the trl definition
      begsample     = event(i).sample*resamp_factor + round(cfg.trialdef.prestim*srate);
      endsample     = event(i).sample*resamp_factor + round(cfg.trialdef.poststim*srate)-1;
      offset        = cfg.trialdef.prestim*srate;  
      trigger       = event(i).value; % remember the trigger (=event type S/F) for each trial
      trl(end+1, :) = [round([begsample endsample offset])  trigger]; 
    end
  end
end