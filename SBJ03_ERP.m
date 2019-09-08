function SBJ03_ERP(SBJ,proc_id,an_id)

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

fig_dir = [root_dir 'PRJ_Error_eeg/results/ERPs/' SBJ '/' an_id '/'];
if ~exist(fig_dir,'dir')
   mkdir(fig_dir);
end

%% Load EEG and Behavioral Data
load([SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat']);

%% Split by condition and outcome
if strcmp(an_id,'dif_out')
ez_idx = strcmp(bhv.Condition,'easy');
hd_idx = strcmp(bhv.Condition,'hard');

cfgs = [];
cfgs.trials = intersect(find(ez_idx),find(bhv.Hit==1));
ez_wn = ft_selectdata(cfgs,trials);
cfgs.trials = intersect(find(ez_idx),find(bhv.Hit==0));
ez_ls = ft_selectdata(cfgs,trials);
cfgs.trials = intersect(find(hd_idx),find(bhv.Hit==1));
hd_wn = ft_selectdata(cfgs,trials);
cfgs.trials = intersect(find(hd_idx),find(bhv.Hit==0));
hd_ls = ft_selectdata(cfgs,trials);

%% ERP Analysis
cfg = [];
ew_erp = ft_timelockanalysis(cfg, ez_wn);
el_erp = ft_timelockanalysis(cfg, ez_ls);
hw_erp = ft_timelockanalysis(cfg, hd_wn);
hl_erp = ft_timelockanalysis(cfg, hd_ls);

%% Plot ERPs
% cfg = [];
% cfg.channels = {'all','-PO7','-O1'};
% cfg.layout    = 'biosemi64.lay';
% cfg.interactive = 'yes';
% cfg.showoutline = 'yes';
% ft_multiplotER(cfg, erp);
for ch_ix = 1:numel(ew_erp.label)
    fig_name = [SBJ '_' an_id '_' ew_erp.label{ch_ix}];
    figure('Name',fig_name);
    subplot(1,2,1); hold on;
    plot(ew_erp.time,ew_erp.avg(ch_ix,:),'k');
    plot(el_erp.time,el_erp.avg(ch_ix,:),'r');
    title('Easy');
    
    subplot(1,2,2); hold on;
    plot(hw_erp.time,hw_erp.avg(ch_ix,:),'k');
    plot(hl_erp.time,hl_erp.avg(ch_ix,:),'r');
    title('Hard');
    
    saveas(gcf,[fig_dir fig_name]);
    close(gcf);
end


%% !!! yet a third script - LOOK AT ERP DIFFERENCE WAVES:s
%do this to not write over orig data?
%subtract the avg of task2 from the average of task 1 at each channel
%cfg = [];
%task1 = ft_timelockanalysis(cfg, comp_clean);
% note that the following appears to do the same:
% difference     = task1;                   % copy one of the structures
% difference.avg = task1.avg - task2.avg;   % compute the difference ERP
% however that will not keep original information, whereas ft_math will

%cfg = [];
%cfg.channels = (1:64);
%cfg.layout    = 'biosemi64.lay';
%cfg.interactive = 'yes';
%cfg.showoutline = 'yes';
%ft_multiplotER(cfg, task1)
%plot the differences at each channel on the specified layout!





