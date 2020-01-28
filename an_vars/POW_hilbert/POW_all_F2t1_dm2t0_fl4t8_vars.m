% Data Selection
an.ROI         = {'all'};             % Channel to be analyzed
an.event_type  = 'F';           % event around which to cut trials
an.trial_lim_s = [-0.25 1];       % window in SEC for cutting trials

% Baseline Correction
an.bsln_lim    = [-0.25 0.05];    % window in SEC for baseline correction
an.bsln_type   = 'zboot';
an.bsln_boots  = 500;

% POW Filtering
an.demean_yn   = 'no';
an.lp_yn       = 'yes';
an.lp_freq     = 8;
an.hp_yn       = 'yes';
an.hp_freq     = 4;
an.hp_filtord  = 1;
an.hilbert     = 'abs';

an.dsamp_yn    = 0;
an.dsamp_freq  = 0;

% Window Logic Check: bsln_lim is within trial_lim_s
if an.bsln_lim(1) < an.trial_lim_s(1) || an.bsln_lim(2) > an.trial_lim_s(2)
    error('an.bsln_lim is outside an.trial_lim_s!');
end

% Wavelet TFR Parameters
% % cfg_tfr.method = 'wavelet'
% % cfg_tfr.output = 'pow';
% % cfg_tfr.taper = 'hanning';
% % cfg_tfr.foi = [2:30]; 
% % cfg_tfr.width = 2; %default
% % cfg_tfr.toi  = 0:0.004:2.8; 

