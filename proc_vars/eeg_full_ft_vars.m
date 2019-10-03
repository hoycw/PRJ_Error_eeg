%% Pipeline Processing Variables: "eeg_full_ft" for EEG analysis
% Trial Cut Parameteres
proc.event_type    = 'S';          % 'S'/'F': lock trial to stim/feedback
proc.trial_lim_s   = [-0.15 2.8];    % data segments (in seconds) to grab around events
if strcmp(proc.event_type,'S')
    proc.event_code = 1;
elseif strcmp(proc.event_type,'F')
    proc.event_code = 2;
else
    error(['Unknown event_type: ' proc.event_type]);
end

% Behavioral Processing
proc.rt_bounds = [0.6 1.4];          % bounds on a reasonable RT to be detected with KLA algorithm
% Varaince-Based Trial Rejection Parameters
proc.var_std_warn_thresh = 3;

% Data Preprocessing
proc.plot_psd      = '1by1';         % type of plot for channel PSDs
proc.resample_yn   = 'yes';
proc.origsample_freq = 1024;
proc.resample_freq = 250;
proc.demean_yn     = 'yes';
proc.reref_yn      = 'yes';
proc.ref_method    = 'avg';
proc.bp_yn         = 'yes';
proc.bp_freq       = [0.1 30];
proc.bp_order      = 2;
proc.hp_yn         = 'no';
proc.hp_freq       = 0.1;            % [] skips this step
proc.hp_order      = 2;              % Leaving blank causes instability error, 1 or 2 works
proc.lp_yn         = 'no';
proc.lp_freq       = 30;             % [] skips this step
proc.notch_yn      = 'no';
proc.notch_type    = 'bandstop';     % method for nothc filtering out line noise

% ICA preprocessing
proc.ICA_hp_yn       = 'no';
proc.ICA_hp_freq     = 0.1;
proc.eog_ic_corr_cut = 0.3;        % EOG IC correlation threshold for tossing ICs % changed to 0.3, ran into a lot of components that were generating errors because not hitting threshold for horizontal eogs
proc.eog_bp_yn       = 'yes';
proc.eog_bp_freq     = [1 15];  % from ft_rejectvisual help page
proc.eog_bp_filtord  = 4;       % from ft_rejectvisual help page
