% Target Time TFR Feature Parameters:
ft.feat_lab    = 'deltaP3';                      % Overall name of features to extract
ft.measure     = 'tfWin';                       
ft.feat_id     = [ft.feat_lab '_' ft.measure];  % unique ID for this feature set
ft.grp_id      = 'DifFB';                        % select group of conditions needed ('rare','Odd','Tar')
ft.an_id       = 'TFR_Pz_F2t1_db2t0_fl1t12';  % TFR analysis ID

ft.name        = {'deltaP3'};                 % name of each feature
% ft.cond not used because features computed for all TT conditions
ft.chan        = {'Pz'};                   % Channel to search for feature
ft.win_lim     = [0.3 0.45];             % Time range to average over
ft.freq_lim    = [1 4];                 % frequency limits (in Hz) to average over


