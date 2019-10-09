plt.plt_lim    = [-0.2 2.8];
plt.x_step_sz  = 0.20;
plt.legend     = 1;
plt.legend_loc_S = 'northwest';
plt.legend_loc_R = 'northeast';
plt.legend_loc = 'northeast';

plt.ylim_fudge = 0.1;

plt.errbar_alpha = 0.5;
%plt.main_style = '--';

%plt.plot_avg  = 0;
%plt.sig_type  = 'patch';
%plt.sig_color = [0.5 0.5 0.5];
%plt.sig_alpha = 0.3;
%plt.sig_style = '-';
%plt.sig_width = 5;
%plt.sig_scat_size  = 70;
%plt.sig_scat_mrkr  = 'o';
%plt.sig_scat_size2 = 30;
%plt.sig_scat_mrkr2 = '+';

plt.evnt_type  = {'S','R','F'};
plt.evnt_width = [2 2 2];
plt.evnt_color = {'k','k','k'};
plt.evnt_style = {'-','--','-'};
plt.clim_perc  = [5, 95];
