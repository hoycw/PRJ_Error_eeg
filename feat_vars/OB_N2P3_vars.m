% ERP Feature Parameters
ft.conditions  = 'rare';                                % select oddball and target trials
ft.feat_name   = {'N2b','P3a','N2c','P3b'};
ft.feat_cond   = {'Odd','Odd','Tar','Tar'};           % trial conditions to model
ft.pk_chan     = {'Fz','Cz','Fz','Pz'};
ft.pk_lim      = [0.2 0.3; 0.3 0.45; 0.2 0.3; 0.3 0.45];
ft.pk_sign     = [-1; 1; -1; 1];
ft.grp_method  = 'sbj';             % {'sbj', 'jackknife'}
% ft.model_id    = [st.model_lab '_' st.model_cond];

