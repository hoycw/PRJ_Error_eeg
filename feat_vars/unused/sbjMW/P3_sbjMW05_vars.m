% Oddball ERP Feature Parameters: N2b and P3a from Oddball condition
ft.feat_lab    = 'P3';                      % Overall name of features to extract
ft.measure     = 'sbjMW';                       % {'sbjPk','sbjMW','grpMW'}
ft.feat_id     = [ft.feat_lab '_' ft.measure];  % unique ID for this feature set
ft.grp_id      = 'DifFB';                        % select group of conditions needed ('rare','Odd','Tar')
ft.an_id       = 'ERP_Pz_F2t1_dm2t0_fl05t20';  % ERP analysis ID

ft.name        = {'P3'};                 % name of each feature
ft.cond        = {'DifFB'};                 % trial condition to search for feature
ft.chan        = {'Pz'};                   % Channel to search for feature
ft.lim         = [0.3 0.45];           % Time range to search for feature
ft.pk_sign     = [1];                       % polarity of the peak
ft.mn_lim      = [-0.025 0.025];      % time limits for mean window (in sec)
ft.plot_lim    = [0 0.6];                       % time limits for plotting ERPs for QA

% ft.grp_method  = 'sbj';             % {'sbj', 'jackknife'}

