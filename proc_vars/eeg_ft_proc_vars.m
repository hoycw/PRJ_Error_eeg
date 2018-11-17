%% Pipeline Processing Variables: "eeg_ft" for EEG analysis

% Data Preprocessing
proc_vars.plot_psd      = '1by1';         % type of plot for channel PSDs
proc_vars.resample_yn   = 'yes';
proc_vars.resample_freq = 250;
proc_vars.demean_yn     = 'yes';
proc_vars.reref_yn      = 'yes';
proc_vars.ref_method    = 'avg';
proc_vars.bp_yn         = 'yes';
proc_vars.bp_freq       = [0.5 40];
proc_vars.hp_yn         = 'no';
%proc_vars.hp_freq       = 0.1;            % [] skips this step
%proc_vars.hp_order      = 4;              % Leaving blank causes instability error, 1 or 2 works
proc_vars.lp_yn         = 'no';
%proc_vars.lp_freq       = 30;            % [] skips this step
proc_vars.notch_type    = 'bandstop';     % method for nothc filtering out line noise

% cleanline   = 'yes';                % Use Tim Mullen's cleanline function
% dft_yn      = 'no';
% bs_yn       = 'no';                % Parameters for this in SBJ_vars

% Behavioral Processing
proc_vars.rt_bounds = [0.3 2.0];          % bounds on a reasonable RT to be detected with KLA algorithm

% Trial Cut Parameteres
proc_vars.event_type    = 'stim';         % 'stim'/'resp': lock trial to these event
proc_vars.trial_lim_s = [-0.25 2];      % data segments (in seconds) to grab around events
proc_vars.RT_std_thresh = 3;              % rejection threshold for RTs

% Varaince-Based Trial Rejection Parameters
proc_vars.var_std_warning_thresh = 3;
