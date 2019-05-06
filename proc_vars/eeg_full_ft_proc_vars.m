%% Pipeline Processing Variables: "eeg_full_ft" for EEG analysis
% Trial Cut Parameteres
proc_vars.event_type    = 'S';          % 'S'/'F': lock trial to stim/feedback
proc_vars.trial_lim_s   = [-0.25 2.8];    % data segments (in seconds) to grab around events
if strcmp(proc_vars.event_type,'S')
    proc_vars.event_code = 1;
elseif strcmp(proc_vars.event_type,'F')
    proc_vars.event_code = 2;
else
    error(['Unknown event_type: ' proc_vars.event_type]);
end

% Behavioral Processing
proc_vars.rt_bounds = [0.6 1.4];          % bounds on a reasonable RT to be detected with KLA algorithm
% Varaince-Based Trial Rejection Parameters
proc_vars.var_std_warn_thresh = 3;

% Data Preprocessing
proc_vars.plot_psd      = '1by1';         % type of plot for channel PSDs
proc_vars.resample_yn   = 'yes';
proc_vars.origsample_freq = 1024;
proc_vars.resample_freq = 250;
proc_vars.demean_yn     = 'yes';
proc_vars.reref_yn      = 'yes';
proc_vars.ref_method    = 'avg';
proc_vars.bp_yn         = 'yes';
proc_vars.bp_freq       = [0.1 30];
proc_vars.bp_order      = 2;
proc_vars.hp_yn         = 'no';
proc_vars.hp_freq       = 0.1;            % [] skips this step
proc_vars.hp_order      = 2;              % Leaving blank causes instability error, 1 or 2 works
proc_vars.lp_yn         = 'no';
proc_vars.lp_freq       = 30;             % [] skips this step
proc_vars.notch_yn      = 'no';
proc_vars.notch_type    = 'bandstop';     % method for nothc filtering out line noise

% ICA preprocessing
proc_vars.ICA_hp_yn       = 'no';
proc_vars.ICA_hp_freq     = 0.1;
proc_vars.eog_ic_corr_cut = 0.3;        % EOG IC correlation threshold for tossing ICs % changed to 0.3, ran into a lot of components that were generating errors because not hitting threshold for horizontal eogs
proc_vars.eog_bp_yn       = 'yes';
proc_vars.eog_bp_freq     = [1 15];  % from ft_rejectvisual help page
proc_vars.eog_bp_filtord  = 4;       % from ft_rejectvisual help page
