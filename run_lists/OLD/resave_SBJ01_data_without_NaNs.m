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
proc_id = 'eeg_full_ft';
SBJs = {'EEG04','EEG06','EEG07','EEG08','EEG10','EEG12','EEG09'};

%% REDO SBJ01 SAVING
for s = 1:numel(SBJs)
    %% Processing variables
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
    eval(proc_vars_cmd);
    
    %% Re-Create data variable
    % Load and preprocess the data
    % Update ear ref channels with prefix/suffix
    if isfield(SBJ_vars.ch_lab,'prefix')
        ear_lab1 = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{1}];
        ear_lab2 = [SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{2}];
    else
        ear_lab1 = [SBJ_vars.ch_lab.ears{1}];
        ear_lab2 = [SBJ_vars.ch_lab.ears{2}];
    end
    if isfield(SBJ_vars.ch_lab,'suffix')
        ear_lab1 = [ear_lab1 SBJ_vars.ch_lab.suffix];
        ear_lab2 = [ear_lab2 SBJ_vars.ch_lab.suffix];
    end
    
    % Load and preprocess
    cfg = [];
    cfg.continuous = 'yes';
    % cfg.lpfilter   = proc.lp_yn;
    % cfg.hpfilter   = proc.hp_yn;
    cfg.bpfilter   = proc.bp_yn;
    cfg.bpfreq     = proc.bp_freq;
    cfg.bpfiltord  = proc.bp_order;
    cfg.demean     = proc.demean_yn;
    cfg.reref      = proc.reref_yn;
    cfg.refmethod  = proc.ref_method;
    cfg.refchannel = {ear_lab1, ear_lab2};
    
    data = cell(size(SBJ_vars.block_name));
    for b_ix = 1:numel(SBJ_vars.block_name)
        cfg.dataset    = SBJ_vars.dirs.raw_filename{b_ix};
        data{b_ix} = ft_preprocessing(cfg);
    end
    data = fn_concat_blocks(data);
    
    % Downsample
    if strcmp(proc.resample_yn,'yes')
        cfg = [];
        cfg.resamplefs = proc.resample_freq;
        data = ft_resampledata(cfg, data);
    end
    
    % Fix channel labels
    % Remove bad channels
    null_neg = cell(size(SBJ_vars.ch_lab.null));
    for null_ix = 1:numel(SBJ_vars.ch_lab.null)
        null_neg{null_ix} = ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.null{null_ix} SBJ_vars.ch_lab.suffix];
    end
    bad_neg = cell(size(SBJ_vars.ch_lab.bad));
    for bad_ix = 1:numel(SBJ_vars.ch_lab.bad)
        bad_neg{bad_ix} = ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.bad{bad_ix} SBJ_vars.ch_lab.suffix];
    end
    rep_neg = cell(size(SBJ_vars.ch_lab.replace));
    for rep_ix = 1:numel(SBJ_vars.ch_lab.replace)
        rep_neg{rep_ix} = ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.replace{rep_ix}{1} SBJ_vars.ch_lab.suffix];
    end
    ears_neg = {['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{1} SBJ_vars.ch_lab.suffix],...
        ['-' SBJ_vars.ch_lab.prefix SBJ_vars.ch_lab.ears{2} SBJ_vars.ch_lab.suffix]};
    
    cfg = [];
    cfg.channel = [{'all'}, bad_neg, rep_neg, ears_neg, null_neg];
    data = ft_selectdata(cfg,data);
    
    % Strip Pre/Suffix if Necessary
    for ch_ix = 1:numel(data.label)
        data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.prefix,'');
        data.label{ch_ix} = strrep(data.label{ch_ix},SBJ_vars.ch_lab.suffix,'');
        % Replace label of the externals with actual electrode name
        for x = 1:numel(SBJ_vars.ch_lab.replace)
            if strcmp(SBJ_vars.ch_lab.replace{x}{2},data.label{ch_ix})
                data.label{ch_ix} = SBJ_vars.ch_lab.replace{x}{1};
            end
        end
    end
    
    % Remove unipolar EOG
    warning('WARNING!!! Assuming Fp2 is the second vertical EOG!');
    eog_v_low_ix = ~strcmp(SBJ_vars.ch_lab.eog_v,'Fp2');    % only toss the lower one
    eog_neg = [fn_ch_lab_negate(SBJ_vars.ch_lab.eog_h),fn_ch_lab_negate(SBJ_vars.ch_lab.eog_v(eog_v_low_ix))];
    trig_neg = fn_ch_lab_negate({SBJ_vars.ch_lab.trigger});
    cfg = [];
    cfg.channel = [{'all'}, eog_neg, trig_neg];
    data = ft_selectdata(cfg,data);
    
    %% Load origianl data
    orig_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
    tmp = load(orig_fname);
    
    % Unpack originals
    eog          = tmp.eog;
    bad_epochs   = tmp.bad_epochs;
    icaunmixing  = tmp.icaunmixing;
    icatopolabel = tmp.icatopolabel;
    
    % Check for differences
    if max(max(data.trial{1}(:, 1:end)-tmp.data.trial{1})) > 0.0001 %%Fixed bug with Size difference!
        error('data is not close enough, double check!');
    end
    
    %% Copy old file
    copy_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '_withNaNs.mat'];
    % check if the copy already exists
    if exist(copy_fname,'file')
        error('copy file already exists, dont overwrite it!');
    end
    copy_cmd = ['cp ' orig_fname ' ' copy_fname];
    system(copy_cmd);
    
    %% Save it out again
    save(orig_fname, 'icaunmixing', 'icatopolabel', 'data', 'eog','bad_epochs');
    
    clear SBJ_vars SBJ icaunmixing icatopolabel data eog bad_epochs
    clear b_ix bad_ix bad_neg cfg ch_ix copy_fname ear_lab1 ear_lab2 ears_neg
    clear eog_neg eog_v_low_ix null_ix null_neg orig_fname proc proc_vars_cmd
    clear rep_ix rep_neg SBJ SBJ_vars_cmd trig_neg x tmp copy_cmd
end
