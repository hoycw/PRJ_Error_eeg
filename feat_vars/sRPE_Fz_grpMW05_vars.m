% Oddball ERP Feature Parameters: N2b and P3a from Oddball condition
ft.feat_lab    = 'sRPE';                      % Overall name of features to extract
ft.measure     = 'grpMW';                       % {'sbjPk','sbjMW','grpMW'}
ft.feat_id     = [ft.feat_lab '_' ft.measure];  % unique ID for this feature set
ft.grp_id      = 'DifFB';                        % select group of conditions needed ('rare','Odd','Tar')
ft.an_id       = 'ERP_Fz_F2t1_dm2t0_fl05t20';  % ERP analysis ID

ft.name        = {'sRPE'};                 % name of each feature
ft.cond        = {'DifFB'};                 % trial condition to search for feature
ft.chan        = {'Fz'};                   % Channel to search for feature
% ft.lim         = [0.3 0.45];           % Time range to search for feature
% ft.pk_sign     = [1];                       % polarity of the peak

ft.pk_reg_id   = 'sRPE';
ft.pk_stat_id  = 'ERPEsL_DifFB_lme_st05t5';            % stat_id from which to get peak time
ft.pk_model_lab = 'ERPEsL';
ft.pk_an_id    = 'ERP_Fz_F2t1_dm2t0_fl05t20';            % an_id from which to get peak time
ft.mn_lim      = [-0.025 0.025];      % time limits for mean window (in sec)
ft.plot_lim    = [0 0.6];                       % time limits for plotting ERPs for QA

% ft.grp_method  = 'sbj';             % {'sbj', 'jackknife'}

