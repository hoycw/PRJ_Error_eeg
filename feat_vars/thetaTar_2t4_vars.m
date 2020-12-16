% Oddball TFR Feature Parameters:
ft.feat_lab    = 'thetaTar';                      % Overall name of features to extract
ft.measure     = 'tfWin';                       
ft.feat_id     = [ft.feat_lab '_' ft.measure];  % unique ID for this feature set
ft.grp_id      = 'Tar';                        % select group of conditions needed ('rare','Odd','Tar')
ft.an_id       = 'TFR_Fz_S2t1_db2t0_fl1t12';  % TFR analysis ID

ft.name        = {'thetaTar'};                 % name of each feature
ft.cond        = {'Tar'}; 
ft.chan        = {'Fz'};                   % Channel to search for feature
ft.win_lim     = [0.2 0.4];             % Time range to average over
ft.freq_lim    = [4 8];                 % frequency limits (in Hz) to average over


