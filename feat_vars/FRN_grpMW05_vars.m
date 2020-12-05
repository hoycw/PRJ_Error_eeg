% Target Time ERP Feature Parameters:
ft.feat_lab    = 'FRN';                      % Overall name of features to extract
ft.measure     = 'grpMW';                       % {'sbjPk','sbjMW','grpMW'}
ft.feat_id     = [ft.feat_lab '_' ft.measure];  % unique ID for this feature set
ft.grp_id      = 'DifFB';                        % select group of conditions needed ('rare','Odd','Tar')
ft.an_id       = 'ERP_Fz_F2t1_dm2t0_fl05t20';  % ERP analysis ID

ft.name        = {'FRN'};                 % name of each feature
% ft.cond not used because features computed for all TT conditions
ft.chan        = {'Fz'};                   % Channel to search for feature
ft.lim         = [0.18 0.3];           % Time range to search for feature
ft.pk_sign     = [-1];                       % polarity of the peak
ft.mn_lim      = [-0.025 0.025];      % time limits for mean window (in sec)
ft.plot_lim    = [0 0.6];                       % time limits for plotting ERPs for QA

% ft.grp_method  = 'sbj';             % {'sbj', 'jackknife'}

