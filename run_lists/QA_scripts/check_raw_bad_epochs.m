SBJ = 'EEG13';
block = 1;
proc_id = 'eeg_full_ft';

%%
root_dir = '/Volumes/hoycw_clust/';
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
load([root_dir 'PRJ_Error_eeg/scripts/utils/cfg_plot_eeg.mat']);

%%
if numel(SBJ_vars.block_name)>1
        block_suffix = ['_' SBJ_vars.block_name{block}];
else
        block_suffix = '';
end
cfg=[];
cfg.dataset  = SBJ_vars.dirs.raw_filename{block};
cfg.demean   = 'yes';
cfg.hpfilter = 'yes';
cfg.hpfreq   = 0.5;
cfg.hpfiltord = 2;
raw = ft_preprocessing(cfg);

% Downsample
if strcmp(proc.resample_yn,'yes')
    cfg = [];
    cfg.resamplefs = proc.resample_freq;
    raw = ft_resampledata(cfg, raw);
end

%%
load([SBJ_vars.dirs.events SBJ '_raw_bad_epochs' block_suffix '.mat']);
cfg_plot.artfctdef.visual.artifact = bad_epochs;

browsed_raw = ft_databrowser(cfg_plot, raw);
