%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

%%
addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% General parameters
SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
           'EEG01','EEG02','EEG03','EEG04','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12'};

%% Linear Mixed Effects Model
proc_id   = 'eeg_full_ft';
an_ids    = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};
stat_id   = 'RL_all_lme_st0t5';
plt_id    = 'ts_F2to1_evnts_sigLine';
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

% for s = 1:numel(SBJs)
%     SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
% end

for an_ix = 1:numel(an_ids)
    SBJ04c_ERP_grp_stats_LME_RL(SBJs,proc_id,an_ids{an_ix},stat_id);
    SBJ04d_ERP_plot_stats_LME_RL(SBJs,proc_id,an_ids{an_ix},stat_id,plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    % SBJ04d_ERP_plot_stats_LME(SBJs,proc_id,an_ids{an_ix},stat_id,plt_id,save_fig,...
    %         'fig_vis',fig_vis,'fig_ftype',fig_ftype,'plot_median',1);
end

%% Compare p values across analyses
proc_id = 'eeg_full_ft';
an_id   = 'ERP_Fz_F2t1_dm2t0_fl05t20';%'ERP_Pz_F2t1_dm2t0_fl05t20';
do_id   = 'DifOut_lme_st0t5';
do_lab  = {'Dif','Out','Dif*Out'};
rl_id   = 'RL_DO_lme_st0t5';
rl_lab  = {'pWin','sPE','uPE'};

% Load DO p values
load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' do_id '_' an_id '.mat']);
do_pvals = nan([numel(do_lab) numel(lme)]);
for grp_ix = 1:numel(do_lab)
    for t_ix = 1:numel(lme)
        do_pvals(grp_ix,t_ix) = lme{t_ix}.Coefficients.pValue(grp_ix+1);
    end
end

% Load RL p values
load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' rl_id '_' an_id '.mat']);
rl_pvals = nan([numel(rl_lab) numel(lme)]);
for grp_ix = 1:numel(rl_lab)
    for t_ix = 1:numel(lme)
        rl_pvals(grp_ix,t_ix) = lme{t_ix}.Coefficients.pValue(grp_ix+1);
    end
end

% Load time vector
eval(['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' do_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m']);
load([SBJ_vars.dirs.proc,SBJs{1},'_',an_id,'.mat']);
cfgs = []; cfgs.latency = st.stat_lim;
roi = ft_selectdata(cfgs, roi);
time_vec = roi.time{1};

% Plot p value comparison
figure;
for grp_ix = 1:3
    subplot(2,3,grp_ix); hold on;
    plot(time_vec,do_pvals(grp_ix,:),'color','r');
    plot(time_vec,rl_pvals(grp_ix,:),'color','k');
    legend(do_lab{grp_ix},rl_lab{grp_ix});
end
for grp_ix = 1:3
    subplot(2,3,3+grp_ix);
    plot(time_vec,do_pvals(grp_ix,:)-rl_pvals(grp_ix,:),'color','b');
    legend('DO - RL');
    ylabel('DO better ... RL better');
end