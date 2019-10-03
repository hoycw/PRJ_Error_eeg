%% Pipeline Processing Variables: "eeg_ft" for EEG analysis
% Trial Cut Parameteres
proc.event_type    = 'F';          % 'S'/'F': lock trial to stim/feedback
proc.trial_lim_s   = [-0.7 1.8];    % data segments (in seconds) to grab around events
if strcmp(proc.event_type,'S')
    proc.event_code = 1;
elseif strcmp(proc.event_type,'F')
    proc.event_code = 2;
else
    error(['Unknown event_type: ' proc.event_type]);
end
%proc.RT_std_thresh = 3;            % rejection threshold for RTs
% Behavioral Processing
proc.rt_bounds = [0.6 1.4];          % bounds on a reasonable RT to be detected with KLA algorithm
% Varaince-Based Trial Rejection Parameters
proc.var_std_warn_thresh = 3;

% Data Preprocessing
proc.plot_psd      = '1by1';         % type of plot for channel PSDs
proc.resample_yn   = 'yes';
proc.resample_freq = 250;
proc.demean_yn     = 'yes';
proc.reref_yn      = 'yes';
proc.ref_method    = 'avg';
proc.bp_yn         = 'yes';
proc.bp_freq       = [0.5 40];
proc.hp_yn         = 'no';
%proc.hp_freq       = 0.1;            % [] skips this step
%proc.hp_order      = 4;              % Leaving blank causes instability error, 1 or 2 works
proc.lp_yn         = 'no';
%proc.lp_freq       = 30;            % [] skips this step
proc.notch_yn      = 'no';
proc.notch_type    = 'bandstop';     % method for nothc filtering out line noise

% cleanline   = 'yes';                % Use Tim Mullen's cleanline function
% dft_yn      = 'no';
% bs_yn       = 'no';                % Parameters for this in SBJ_vars

proc.eog_ic_corr_cut = 0.3;        % EOG IC correlation threshold for tossing ICs
