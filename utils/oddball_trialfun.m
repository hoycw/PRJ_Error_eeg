function trl = oddball_trialfun(cfg)
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
if cfg.blocknum > 1
   index_overall = cfg.endb1;
else
    index_overall = 1;
end
%checks to see if starting with oddball or not, oddball_section = 1 means yes
oddball_section = 0;
if index_overall < 400
    for i = index_overall:400
        if (event(i).value == 3)
            oddball_section = 1;
        end
    end
end
% Find event cuts and built trl matrix
trl = [];
for i=1:length(event)
  if strcmp(event(i).type, cfg.trialdef.eventtype)
    if (event(i).value == 255 || event(i).value == 254) && i == max(cfg.tt_trigger_ix, cfg.odd_trigger_ix) %the second paradigm's index
        % 255 and 254 seem to be interchangable.
        % This marks the end of the oddball section and the start of the TT
        oddball_section = abs(oddball_section - 1); % I think this is a cheat way to swithc between 0 and 1
        if oddball_section == 1
            oddball_sample = event(i).sample;
            tt_sample = 1;
        end
        if oddball_section == 0
            tt_sample = event(i).sample;
            oddball_sample = 1;
        end
    end
    if (event(i).value == 255 || event(i).value == 254) && oddball_section
        % This restarts the trl if you had to restart the paradigm and lsot
        % the response log
        trl = [];
    end
     % Add trials in the oddball section with real event codes (1,2,3)
     % Add oddball stim onsets: real event codes (1,2,3) = (std,tar,odd)
     if oddball_section && (event(i).value~=254 && event(i).value~=255)
       begsample     = event(i).sample*resamp_factor + round(cfg.trialdef.prestim*srate);
       endsample     = event(i).sample*resamp_factor + round(cfg.trialdef.poststim*srate)-1;
      offset        = cfg.trialdef.prestim*srate;  
      trigger       = event(i).value; % remember the trigger (=condition) for each trial
      trl(end+1, :) = [round([begsample endsample offset])  trigger]; 
    end
  end
end
