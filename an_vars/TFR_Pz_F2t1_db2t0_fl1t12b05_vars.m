% Data Selection
an.ROI         = {'Pz'};             % Channel to be analyzed
an.event_type  = 'F';           % event around which to cut trials
an.trial_lim_s = [-0.2 1];       % window in SEC for cutting trials
an.bsln_type   = 'db';
an.bsln_lim    = [-0.2 0];
an.bsln_boots  = 0;
% an.demean_yn   = 'no';    % don't need this for TFRs?
an.avgoverfreq = 0;
an.itpc        = 0;

% TFR Parameters
cfg_tfr = [];
cfg_tfr.method     = 'wavelet';
cfg_tfr.output     = 'pow';
cfg_tfr.taper      = 'hanning';
cfg_tfr.foi        = [1:0.5:4 5:12];
cfg_tfr.width      = 3; %default
cfg_tfr.toi        = 'all'; %-0.2:0.004:1.0;
cfg_tfr.keeptrials = 'yes'; % need trials for stats, can average later
cfg_tfr.keeptapers = 'no';
cfg_tfr.pad        = 'maxperlen';                         %add time on either side of window
cfg_tfr.padtype    = 'zero';

% Window Logic Check: bsln_lim is within trial_lim_s
if an.bsln_lim(1) < an.trial_lim_s(1) || an.bsln_lim(2) > an.trial_lim_s(2)
    error('an.bsln_lim is outside an.trial_lim_s!');
end
