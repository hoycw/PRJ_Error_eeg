function trl = tt_trialfun(cfg);

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim
% cfg.dsrate (sampling rate for downsampling, just write orig (1024) if no
% downsampling

hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);
fprintf('%i events found!\n',numel(event));
%Changed this to allow for downsampling
trl = [];
 
for i=1:length(event)
  if strcmp(event(i).type, cfg.trialdef.eventtype)
    % it is a trigger, see whether it has the right value
    if ismember(event(i).value, cfg.trialdef.eventvalue)
      % add this to the trl definition
      begsample     = event(i).sample/hdr.Fs*cfg.dsrate + round(cfg.trialdef.prestim*cfg.dsrate);
      endsample     = event(i).sample/hdr.Fs*cfg.dsrate + round(cfg.trialdef.poststim*cfg.dsrate)-1;
      offset        = cfg.trialdef.prestim*cfg.dsrate;  
      trigger       = event(i).value; % remember the trigger (=condition) for each trial
      trl(end+1, :) = [round([begsample endsample offset])  trigger]; 
    end
  end
end