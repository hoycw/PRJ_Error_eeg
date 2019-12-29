function output = tt_trialsamples(cfg);
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
output.orig_sample = hdr.Fs;
%checks to see if starting with oddball or not, oddball_section = 1 means yes
oddball_section = 0;
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
    if (event(i).value == 255 || event(i).value == 254) && i == max(cfg.tt_trigger_ix, cfg.odd_trigger_ix)
        % 255 and 254 seem to be interchangable.
        oddball_section = abs(oddball_section - 1); %I think this is a cheat way to flip between 0 and 1
        if oddball_section == 1
            output.oddball_sample = event(i).sample;
            output.tt_sample = 1;
        end
        if oddball_section == 0
            output.tt_sample = event(i).sample;
            output.oddball_sample = 1;
        end
    end
  end
end

end