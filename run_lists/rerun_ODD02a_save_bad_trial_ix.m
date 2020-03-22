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
proc_id = 'odd_full_ft';
SBJ_id = 'goodEEG';

sbj_file = fopen([root_dir 'PRJ_Error_EEG/scripts/SBJ_lists/' SBJ_id '.sbj']);
tmp = textscan(sbj_file,'%s');
fclose(sbj_file);
SBJs = tmp{1}; clear tmp;

proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

%% REDO SBJ01 SAVING
for s = 1:numel(SBJs)
    %% Processing variables
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% Load data
%     if numel(SBJ_vars.block_name)>2; error('not ready for 3+ blocks!'); end
    
    % Load EEG
    load([SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'],'data','bad_epochs');
    
    % Load Behavior
    [bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);
    
    %% Cut into trials
    % Need to recut trials on updated data with the nans
    cfg_trl_unconcat = cell(size(SBJ_vars.block_name));
    for b_ix = 1:numel(SBJ_vars.block_name)
        cfg = [];
        cfg.dataset             = SBJ_vars.dirs.raw_filename{b_ix};
        cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
        cfg.trialdef.eventvalue = proc.event_code;        % feedback cocde
        cfg.trialdef.prestim    = proc.trial_lim_s(1);
        cfg.trialdef.poststim   = proc.trial_lim_s(2);
        if startsWith(SBJ, 'EEG')
            cfg.tt_trigger_ix       = SBJ_vars.tt_trigger_ix;
            cfg.odd_trigger_ix      = SBJ_vars.odd_trigger_ix;
        end
        cfg.trialfun            = 'tt_trialfun';
        % Add downsample frequency since triggers are loaded from raw file
        cfg.resamp_freq         = proc.resample_freq;
        cfg_trl_unconcat{b_ix}  = ft_definetrial(cfg);
    end
    
    % Concatenate event times across blocks after adjusting times for gaps
    cfg_trl = cfg_trl_unconcat{1};
    if numel(SBJ_vars.block_name)>1
        % Get length of first block
        hdr = ft_read_header(SBJ_vars.dirs.raw_filename{1});
        endsample = hdr.nSamples;
        origFs = hdr.Fs;
        
        % Adjust event times for concatenated blocks
        cfg_trl_unconcat{2}.trl(:,1) = cfg_trl_unconcat{2}.trl(:,1)+endsample/origFs*proc.resample_freq;
        cfg_trl_unconcat{2}.trl(:,2) = cfg_trl_unconcat{2}.trl(:,2)+endsample/origFs*proc.resample_freq;
        cfg_trl.trl = vertcat(cfg_trl.trl, cfg_trl_unconcat{2}.trl);
    end
    
    % If the recording was started part way through, toss events not recorded
    if any(cfg_trl.trl(:,1)<1)
        cfg_trl.trl(cfg_trl.trl(:,1)<1,:) = [];
    end
    event_onsets = cfg_trl.trl(:,1)-cfg_trl.trl(:,3);
    
    % Cut the data into trials
    trials = ft_redefinetrial(cfg_trl,data);
    
    % Check that behavioral and EEG event triggers line up
    if (numel(bhv.trl_n))~=numel(event_onsets)
        error(['Mismatch in behavioral and neural trial counts: ' num2str((numel(bhv.trl_n)))...
            ' behavioral; ' num2str(numel(event_onsets)) ' neural']);
    end
    
    %% Exclude bad_trials
    % Find trials that overlap with bad_epochs from raw visual inspection
    if ~isempty(bad_epochs)
        bad_raw_trials = fn_find_trials_overlap_epochs(bad_epochs,1:size(data.trial{1},2),...
            event_onsets,proc.trial_lim_s*data.fsample);
    else
        bad_raw_trials = [];
    end
    
    % Identify training and bad behavioral trials
    training_ix = find(bhv.blk==0);
    rt_low_ix   = find(bhv.rt <= proc.rt_bounds(1));
    rt_high_ix  = find(bhv.rt >= proc.rt_bounds(2));
    exclude_trials = unique(vertcat(bad_raw_trials, training_ix, rt_low_ix, rt_high_ix));
    fprintf(2,'\tWarning: Removing %d trials (%d bad_raw, %d training, %d rts)\n', numel(exclude_trials),...
        numel(bad_raw_trials), numel(training_ix), numel(rt_low_ix)+numel(rt_high_ix));
    
    %% Save excluded trial indices
    out_fname = [SBJ_vars.dirs.events SBJ '_' proc_id '_02a_orig_exclude_trial_ix.mat'];
    save(out_fname,'-v7.3','bad_raw_trials','training_ix','rt_low_ix','rt_high_ix','exclude_trials');
    
    clear SBJ_vars SBJ SBJ_vars_cmd data bhv bad_epochs cfg cfg_trl cfg_trl_unconcat
    clear cfg cfg_trl cfg_trl_unconcat event_onsets trials b_ix out_fname
    clear bad_raw_trials training_ix rt_low_ix rt_high_ix exclude_trials
end
