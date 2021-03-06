function SBJ02a_artifact_rejection(SBJ, proc_id, gen_figs, fig_vis, clear_plot)
%% Reject EOG artifacts and bad behavior
% (1) Segment trials
% (2) Cut out raw bad_epoch, training, and bad RT trials
% (3) Identify EOG ICs via correlation
% (4) Generate figures for quality checks on ERPs and ICA components.
%   QA figures are single trial stacks + mean (across all conditions) for both evoked and ICA data
% INPUTS:
%   SBJ [str] - name of the SBJ
%   proc_id [str] - name of the preprocessing pipeline parameters (e.g., 'egg_full_ft')
%   gen_figs [0/1] - binary flag to generate QA figures or not
%   fig_vis [0/1] - if a data_browser view of the time course of the ICA
%   clear_plot [0/1] - binary flag to delete previous QA plot directory components if desired
% OUTPUTS: EEG/ICA/EOG
%   trials [FT struct] - EEG data cut to good trials
%   cfg_trl [struct] - FT trl structure used to cut trials from full data
%   ica [FT struct] - reconstructed ICA data
%   heog_ics [array] - indices of ICs matched to horizontal EOG
%   veog_ics [array] - indices of ICs matched to vertical EOG
%   eog_trials [FT struct] - EOG data cut to good trials
% OUTPUTS: BEHAVIOR
%   bhv [struct] - behavioral data excluding bad trials
%   bhv_fields [cell array] - list of fields in bhv struct
% OUTPUTS: TRIAL EXCLUSION TRACKING
%   NOTE: saving original trial indices means the correct trials can be tossed when
%       loading raw data, e.g., in SBJ05 for TFR filtering on whole data
%   bad_raw_trials [array] - original trial indices tossed for visual artifact marking in SBJ00
%   training_ix [array] - original trial indices tossed for training
%   rt_low_ix [array] - original trial indices tossed for too short RTs
%   rt_high_ix [array] - original trial indices tossed for too long RTs
%   exclude_trials [array] - original trial indices tossed for any reason (all unique trial indices combined)

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults
% NOTE: Colin deleted below line (and Sheila created directory/functions) on 8/14/20,
%   hopefully because it's not necessary anymore (didn't run to check)
% addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);

%% Processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

%% Load data
if numel(SBJ_vars.block_name)>2; error('not ready for 3+ blocks!'); end

% Load EEG
load([SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat']);

% Load Behavior
[bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);

%% Cut into trials
% Need to recut trials on updated data with the nans
cfg_trl_unconcat = cell(size(SBJ_vars.block_name));
for b_ix = 1:numel(SBJ_vars.block_name)
    cfg = [];
    cfg.dataset             = SBJ_vars.dirs.raw_filename{b_ix};
    cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;    % Name of channel with triggers
    cfg.trialdef.eventvalue = proc.event_code;          % feedback cocde
    cfg.trialdef.prestim    = proc.trial_lim_s(1);      % pre-trigger epoch limit
    cfg.trialdef.poststim   = proc.trial_lim_s(2);      % post-trigger epoch limit
    cfg.tt_trigger_ix       = SBJ_vars.tt_trigger_ix;   % trigger code for start of Target Time task
    cfg.odd_trigger_ix      = SBJ_vars.odd_trigger_ix;  % trigger code for start of Oddball task
    if b_ix > 1
        % Check above limits this to only b_ix == 2
        cfg.endb1 = cfg_trl_unconcat{b_ix - 1}.endb1;   % end sample of first data recording
    end
    % NOTE: Colin commented the 4 lines below on 8/14/20 in code review
    %   Not necessary due to above lines, presumably this field was added
    %   to EP* SBJs at some point?
%     if startsWith(SBJ, 'EEG')
%         cfg.tt_trigger_ix       = SBJ_vars.tt_trigger_ix;
%         cfg.odd_trigger_ix      = SBJ_vars.odd_trigger_ix;
%     end
    cfg.trialfun            = 'tt_trialfun';            % custom trial segmentation function
    % Add downsample frequency since triggers are loaded from raw file
    cfg.resamp_freq         = proc.resample_freq;
    cfg.blocknum            = b_ix;
    cfg_trl_unconcat{b_ix}  = ft_definetrial(cfg);
    if b_ix == 1
        cfg_trl_unconcat{b_ix}.endb1 = length(cfg_trl_unconcat{b_ix}.trl);
    end
end

% Concatenate trl struct event times across blocks after adjusting times for gaps
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
eog_trials = ft_redefinetrial(cfg_trl,eog);

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
training_ix = find(bhv.blk==0);                 % Training block index is 0 here, -1 in python
rt_low_ix   = find(bhv.rt <= proc.rt_bounds(1));
rt_high_ix  = find(bhv.rt >= proc.rt_bounds(2));
exclude_trials = unique(vertcat(bad_raw_trials, training_ix, rt_low_ix, rt_high_ix));
fprintf(2,'\tWarning: Removing %i trials (%i bad_raw, %i training, %i rts)\n', numel(exclude_trials),...
    numel(bad_raw_trials), numel(training_ix), numel(rt_low_ix)+numel(rt_high_ix));

% Exclude bad trials
cfgs = [];
cfgs.trials = setdiff([1:numel(trials.trial)], exclude_trials');
trials = ft_selectdata(cfgs, trials);
eog_trials = ft_selectdata(cfgs, eog_trials);

bhv_fields = fieldnames(bhv);
for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(exclude_trials) = [];
end

%% EOG vs. ICA Correlation
% Rebuild the components from trial cut data
cfg           = [];
cfg.unmixing  = icaunmixing;
cfg.topolabel = icatopolabel;
ica           = ft_componentanalysis(cfg, trials);

% Filter EOG
if strcmp(proc.eog_bp_yn,'yes')
    cfg           = [];
    cfg.bpfilter  = proc.eog_bp_yn;
    cfg.bpfreq    = proc.eog_bp_freq;
    cfg.bpfiltord = proc.eog_bp_filtord;
    eog = ft_preprocessing(cfg,eog_trials);
end

% Correlate EOG with ICs
eog_ic_corr = zeros([numel(eog.label), numel(ica.topolabel), numel(eog.trial)]);
for eog_ix = 1:2
    for ic_ix = 1:numel(ica.label)
        for t_ix = 1:numel(ica.trial)
            temp = corrcoef(eog.trial{t_ix}(eog_ix,:), ica.trial{t_ix}(ic_ix,:));
            eog_ic_corr(eog_ix, ic_ix, t_ix) = temp(1,2);
        end
    end
end

% Identify components correlated with EOG above threshold
avg_eog_ic_corr = mean(eog_ic_corr,3);
heog_ics = find(abs(avg_eog_ic_corr(1,:))>proc.eog_ic_corr_cut);
veog_ics = find(abs(avg_eog_ic_corr(2,:))>proc.eog_ic_corr_cut);
if all([isempty(heog_ics), isempty(veog_ics)])
    error('No EOG ICs found!');
elseif isempty(heog_ics); warning('No HEOG IC found!');
elseif isempty(veog_ics); warning('No VEOG IC found!');
end

%% Generate Figures
if gen_figs
    if clear_plot
       dir_name = [SBJ_vars.dirs.proc 'plot/*'];
       delete(dir_name);
    end
    % Plot EOG-ICA Correlations
    figure('Visible',1); hold on;
    scatter(avg_eog_ic_corr(1,:),avg_eog_ic_corr(2,:));
    scatter(avg_eog_ic_corr(1,heog_ics),avg_eog_ic_corr(2,heog_ics),'filled','r');
    scatter(avg_eog_ic_corr(1,veog_ics),avg_eog_ic_corr(2,veog_ics),'filled','r');
    title('EOG ICA Correlations');
    xlabel('HEOG (r)');
    ylabel('VEOG (r)');
    legend('non-sig','sig','location','best');
    eog_ic_fname = [SBJ_vars.dirs.preproc SBJ '_eog_ic_corr_' proc_id '.png'];
    saveas(gcf,eog_ic_fname);
    
    % Print EOG-ICA correlation results
    fprintf('==================== EOG vs. ICA Correlations =========================\n');
    fprintf(['%i ICs correlate with HEOG, r = ' repmat('%.02f ',[1 numel(heog_ics)]) '\n'],...
        numel(heog_ics),avg_eog_ic_corr(1,heog_ics));
    fprintf(['%i ICs correlate with VEOG, r = ' repmat('%.02f ',[1 numel(veog_ics)]) '\n'],...
        numel(veog_ics),avg_eog_ic_corr(2,veog_ics));
    fprintf('=======================================================================\n');
    
    % View ICA Power Spectrum, Variance, Topo
    cfg = [];
    cfg.layout   = 'biosemi64.lay';
    cfg.channel  = 'all';
    cfg.path     = SBJ_vars.dirs.proc;
    cfg.prefix   = 'ICA';
    cfg.viewmode = 'component';
    cfg.fig_vis  = fig_vis;
    fn_icabrowser_modified(SBJ, cfg, ica);
    
    % Plot IC single trial stacks + ERPs
    fn_plot_ERP_stack(SBJ, proc_id, 'ERPstack_full_evnts', ica, 'off', 1, cfg.path);
end
% Plot IC in ft_databrowser
if strcmp(fig_vis,'on')
    cfg = [];
    cfg.viewmode = 'component';
    cfg.channel = 'all';
    cfg.layout   = 'biosemi64.lay';
    ft_databrowser(cfg, ica);
end

%% Save Data
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat'];
save(clean_data_fname, '-v7.3', 'trials', 'cfg_trl', 'ica', 'heog_ics', 'veog_ics', 'eog_trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_02a.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv', 'bhv_fields');

% Save excluded trial indices
excluded_fname = [SBJ_vars.dirs.events SBJ '_' proc_id '_02a_orig_exclude_trial_ix.mat'];
save(excluded_fname,'-v7.3','bad_raw_trials','training_ix','rt_low_ix','rt_high_ix','exclude_trials');

end
