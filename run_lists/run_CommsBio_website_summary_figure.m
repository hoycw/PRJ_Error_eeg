%% Summary figure for Sequential PE Manuscript at Communications Biology
%   This run script creates the summary figure for the Communications
%   Biology website showing the model coefficient topographies for:
%       1) FRN: RPE magnitude in negative valence conditions at 216 ms
%       2) RewP: RPE magnitude in positive valence conditions at 308 ms
%       3) Late Positivity: Probability in negative valence conditions at 380 ms
%           NOTE: Negative chosen since the topo matches the pattern for
%           all conditions better than for positive valence
% Code written by Colin W. Hoy

%% Beta Topographies: Linear Mixed Effects Model (Mean Windows)
SBJ_id    = 'goodall';
an_id     = 'ERP_all_F2t1_dm2t0_fl05t20';
stat_ids  = {'uRPEL_Neg_lme_mn05man216','uRPEL_Pos_lme_mn05man308','uRPEL_Neg_lme_mn05man380'};
plt_id    = 'topo_F18t25';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'svg';

% Plot Beta Topo across time points
SBJ04d_ERP_plot_stats_LME_RL_topo_ts_reg(SBJ_id,an_id,stat_ids,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype);

