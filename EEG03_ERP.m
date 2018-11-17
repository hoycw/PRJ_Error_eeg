function erpscript(filename, SBJ)
addpath(genpath('/Users/SCS22/Desktop/Knight_Lab/Preprocessing_Work/'));
data_in_filename = [SBJ, 'cleaned']
load(data_in_filename, 'data');
SBJ_vars_cmd = ['run ' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
%data_dir = '/Users/SCS22/Desktop/Pilot2_ICATesting/pilot01/';
%data_filename = [data_dir 'savedpilot01data'];
%data_out_filename = [data_dir 'erp_pilot01'];
%load(data_filename, 'data');
%% erp analysis
cfg = [];
erps = ft_timelockanalysis(cfg, data);
cfg = [];
cfg.channels = {'all', '-EXG*'};
cfg.layout    = 'biosemi64.lay';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, erps)
savefig(data_out_filename_erp);
save(data_out_filename_erp, 'erps');