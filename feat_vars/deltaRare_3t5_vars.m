% Oddball TFR Feature Parameters:
ft.feat_lab    = 'deltaRare';                      % Overall name of features to extract
ft.measure     = 'tfWin';                       
ft.feat_id     = [ft.feat_lab '_' ft.measure];  % unique ID for this feature set
ft.grp_id      = 'rare';                        % select group of conditions needed ('rare','Odd','Tar')
ft.an_id       = 'TFR_Pz_S2t1_db2t0_fl1t12';  % TFR analysis ID

ft.name        = {'deltaOdd','deltaTar'};                 % name of each feature
ft.cond        = {'Odd','Tar'}; 
ft.chan        = {'Pz','Pz'};                   % Channel to search for feature
ft.win_lim     = [0.3 0.5; 0.3 0.5];             % Time range to average over
ft.freq_lim    = [1 4; 1 4];                 % frequency limits (in Hz) to average over


