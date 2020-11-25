% ERP Feature Parameters
ft.feat_name   = {'N2b','P3a','N2c','P3b'};
ft.feat_cond   = {'Odd','Odd','Tar','Tar'};           % trial conditions to model
ft.model_id    = [st.model_lab '_' st.model_cond];

% Peak selection parameters
ft.pk_lim      = [0.1 0.26; 0.18 0.3];
ft.pk_sign     = [1; -1];
ft.grp_method  = 'sbj';             % {'sbj', 'jackknife'}

