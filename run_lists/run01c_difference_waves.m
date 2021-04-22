%% ERP difference wave computations and plotting for Sequential PE Revision 1 Submission
% Developed over time, but editted 4/21/21 by Colin W Hoy
%   Sup. Fig. 2: SBJ03c_ERP_plot_diff_grp

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
SBJ_id = 'goodall';%'good1';%
SBJs = fn_load_SBJ_list(SBJ_id);

%% ERP Difference Waves: Fz and Pz over time
%   These condition contrasts are designed to compare multiple regresion results to classic difference waves
an_ids     = {'ERP_Fz_F2t1_dm2t0_fl05t20','ERP_Pz_F2t1_dm2t0_fl05t20'};%
conditions = {'RewP','Pos-Neg','Large-Small','Unlik-Lik'};
proc_id    = 'eeg_full_ft';
save_fig   = 1;
fig_vis    = 'on';
fig_ftype  = 'png';%

for an_ix = 1:numel(an_ids)
    for diff_ix = 1:numel(conditions)
%         plt_id     = 'stack_F2t1_evnt_c5';
        for s = 1:numel(SBJs)
            % Plot SBJ difference waves
%             SBJ03b_ERP_plot_diff(SBJs{s},conditions{diff_ix},proc_id,an_ids{an_ix},plt_id,save_fig,...
%                 'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        end
        
        % Plot Group difference waves
        %*** plots Sup. Fig. 2
        plt_id = 'ts_F2t8_evnts_sigLine';
        SBJ03c_ERP_plot_diff_grp(SBJ_id,conditions{diff_ix},proc_id,an_ids{an_ix},plt_id,save_fig,...
            'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        
        % Plot all SBJ ERPs overlapping (butterfly)
        %     plt_id = 'ts_F2to1_but_evnts_sigPatch';
        %     SBJ03c_ERP_plot_grp_butterfly(SBJs,conditions,proc_id,an_ids{an_ix},plt_id,save_fig,...
        %         'fig_vis',fig_vis,'fig_ftype',fig_ftype);
        %     close all;
    end
end

