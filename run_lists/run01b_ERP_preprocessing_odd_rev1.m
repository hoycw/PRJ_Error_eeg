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
SBJ_id = 'good2';%'goodEEG1';
SBJs = fn_load_SBJ_list(SBJ_id);

%% Run preprocessing
% proc_id = 'odd_full_ft';
% proc_id_ica = 'odd_full_ft';
% gen_figs    = 0;
% fig_vis     = 'off';
% reject_visual = 0;
% plot_final_check = 0;
% 
% SBJ_times = zeros(size(SBJs));
% tic;
% for s = 1:numel(SBJs)
%     SBJ00_raw_view(SBJs{s}, 1, proc_id, 1)
%     ODD01_preproc(SBJs{s}, proc_id)
%     ODD02a_artifact_rejection(SBJs{s}, proc_id, gen_figs, fig_vis)
%     SBJ02b_ica_rejection(SBJs{s}, proc_id, proc_id_ica, reject_visual);
%     SBJ02c_trial_rejection(SBJs{s}, proc_id, plot_final_check)
%     SBJ_times(s) = toc;
%     if s==1; elapsed = SBJ_times(s); else; elapsed = SBJ_times(s)-SBJ_times(s-1); end
%     fprintf('%s preprocessed at %.1f s (SBJ time = %.1f)\n',SBJs{s},SBJ_times(s),elapsed);
% end

%% ERPs: Fz and Pz
proc_id    = 'odd_full_ft';
an_ids     = {'ERP_Fz_S2t1_dm2t0_fl05t20','ERP_Pz_S2t1_dm2t0_fl05t20'};
conditions = 'OB';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for an_ix = 1:numel(an_ids)
    for s = 1:numel(SBJs)
%         SBJ03a_ERP_save(SBJs{s},proc_id,an_ids{an_ix});
%         plt_id     = 'ts_S2t1_evnts_sigLine';
%         SBJ03b_ERP_plot(SBJs{s},conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
%             'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    end
    
    % Group ERP plot
    plt_id     = 'ts_S2t1_evnts_sigLine';
    SBJ03c_ERP_plot_grp(SBJ_id,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
        'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    % plt_id = 'ts_F2to13_but_evnts_sigPatch';
    % SBJ03c_ERP_plot_grp_butterfly(SBJs,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
    %     'fig_vis',fig_vis,'fig_ftype',fig_ftype);
    
    close all;
end

%% ERPs: Full Cap Topography
proc_id    = 'odd_full_ft';
an_id      = 'ERP_all_S2t1_dm2t0_fl05t20';
conditions = 'OB';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';

for s = 1:numel(SBJs)
%     SBJ03a_ERP_save(SBJs{s},proc_id,an_id);
end

% Group Topo Plot: P3
plt_id = 'topo_F3t45';
SBJ03c_ERP_plot_grp_topo_cond(SBJ_id,conditions,proc_id,an_id,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype);

% Group Topo Plot: N2
plt_id = 'topo_F2t3';
SBJ03c_ERP_plot_grp_topo_cond(SBJ_id,conditions,proc_id,an_id,plt_id,save_fig,...
    'fig_vis',fig_vis,'fig_ftype',fig_ftype);

%% Save topo for CPA prototype selection
proc_id    = 'odd_full_ft';
an_id      = 'ERP_all_S2t1_dm2t0_fl05t20';
conditions = 'OB';

SBJ03c_ERP_save_grp_topo_cond(SBJ_id,conditions,proc_id,an_id);

