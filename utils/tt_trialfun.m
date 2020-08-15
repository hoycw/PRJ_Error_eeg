function trl = tt_trialfun(cfg)
% Trial cutting function for Target Time EEG experiment
%   deals with possible Oddball task run before Target Time task
%   Oddball should be run first, but was run second on at least one dataset (EEG07)
%   COLIN 8/14/20: I don't totally understand the task switch detection
%       logic without looking into it more...
% Required Fields:
%   cfg.dataset [str] - pathname to dataset from which to read the events
%   cfg.trialdef.eventtype [str] - name of the channel with the event codes
%   cfg.trialdef.eventvalue [int] - code for desired segmenting events
%   cfg.trialdef.prestim [float] - latency in sec to cut pre-event
%   cfg.trialdef.poststim [float] - latency in sec to cut post-event
% Optional Fields:
%   cfg.resamp_freq [int] - sampling rate if data was resampled (e.g., downsampling)

% Read header and event triggers
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset, 'header', hdr);
fprintf('%i events found!\n',numel(event));

% Set starting index
if cfg.blocknum > 1
    index_overall = cfg.endb1;
else
    index_overall = 1;
end

% Check to see if starting with oddball or not, oddball_section = 1 means yes
oddball_section = 0;
if index_overall < 400
    % Scan early events for 3rd trigger type, which indicates oddball task
    % There should be ~400 events for oddball task
    for i = index_overall:400
        if (event(i).value == 3)
            oddball_section = 1;
        end
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

% Find event cuts and build trl matrix
trl = [];
for i=1:length(event)
    if strcmp(event(i).type, cfg.trialdef.eventtype)
        % Detect task switch, determine which task came second
        if (event(i).value == 255 || event(i).value == 254) && i == max(cfg.tt_trigger_ix, cfg.odd_trigger_ix)
            % Sheila original comment: 255 and 254 seem to be interchangable.
            % COLIN 8/14/20: 254 should be start, and 255 should be end, but unclear...
            oddball_section = abs(oddball_section - 1); %I think this is a cheat way to flip between 0 and 1
            if oddball_section == 1
                cfg.oddball_sample = event(i).sample;
                cfg.tt_sample = 1;
            end
            if oddball_section == 0
                cfg.tt_sample = event(i).sample;
                cfg.oddball_sample = 1;
            end
        end
        if event(i).value == 254
            %if event(i).value == 254 || event(i).value == 255 && ~oddball_section
            trl = [];
        end
        % Add TT events if event code matches eventvalue: 1 (stim) or 2 (feedback)
        if ~oddball_section && ismember(event(i).value, cfg.trialdef.eventvalue)
            % add this to the trl definition
            begsample     = event(i).sample*resamp_factor + round(cfg.trialdef.prestim*srate);
            endsample     = event(i).sample*resamp_factor + round(cfg.trialdef.poststim*srate)-1;
            offset        = cfg.trialdef.prestim*srate;
            trigger       = event(i).value; % remember the trigger type (stim or feedback) for each trial
            trl(end+1, :) = [round([begsample endsample offset])  trigger];
        end
    end
end

end