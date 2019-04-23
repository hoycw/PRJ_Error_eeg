function EEG02a_artifact_rejection(SBJ, proc_id, gen_figs, fig_vis)

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);

%% Load data
% Load EEG
data_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
load(data_fname);

% Load Behavior
[bhv, bhv_fields] = fn_read_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);

%% Cut into trials
% Need to recut trials on updated data with the nans
cfg = [];
cfg.dataset             = SBJ_vars.dirs.raw_filename;
cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
cfg.trialdef.eventvalue = proc_vars.event_code;        % feedback cocde
cfg.trialdef.prestim    = proc_vars.trial_lim_s(1);
cfg.trialdef.poststim   = proc_vars.trial_lim_s(2);
cfg.trialfun            = 'tt_trialfun';
% Add downsample frequency since triggers are loaded from raw file
cfg.resamp_freq         = proc_vars.resample_freq;
cfg_trl = ft_definetrial(cfg);

% If the recording was started part way through, toss events not recorded
if any(cfg_trl.trl(:,1)<1)
    cfg_trl.trl(cfg_trl.trl(:,1)<1,:) = [];
end
event_onsets = cfg_trl.trl(:,1)-cfg_trl.trl(:,3);


% Check that behavioral and EEG event triggers line up
if numel(bhv.Trial)~=numel(event_onsets)
    error(['Mismatch in behavioral and neural trial counts: ' num2str(numel(bhv.Trial))...
        ' behavioral; ' num2str(numel(event_onsets)) ' neural']);
end

% Cut the data into trials
trials = ft_redefinetrial(cfg_trl,data);
eog_trials = ft_redefinetrial(cfg_trl,eog);

%% Exclude bad_trials
% Find trials that overlap with bad_epochs from raw visual inspection
load([SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat']);
if ~isempty(bad_epochs)
    bad_raw_trials = fn_find_trials_overlap_epochs(bad_epochs,1:size(data.trial{1},2),...
        event_onsets,proc_vars.trial_lim_s*data.fsample);    
else
    bad_raw_trials = [];
end

% Identify training and bad behavioral trials
training_ix = find(bhv.Block==-1);
rt_low_ix   = find(bhv.RT <= proc_vars.rt_bounds(1));
rt_high_ix  = find(bhv.RT >= proc_vars.rt_bounds(2));
exclude_trials = unique(vertcat(bad_raw_trials, training_ix, rt_low_ix, rt_high_ix));

% Exclude bad trials
cfgs = [];
cfgs.trials = setdiff([1:numel(trials.trial)], exclude_trials');
trials = ft_selectdata(cfgs, trials);
eog_trials = ft_selectdata(cfgs, eog_trials);

for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(exclude_trials) = [];
end

%% EOG vs. ICA Correlation
% Rebuild the components
cfg           = [];
cfg.unmixing  = icaunmixing;
cfg.topolabel = icatopolabel;
ica           = ft_componentanalysis(cfg, trials);

% Filter EOG
if strcmp(proc_vars.eog_bp_yn,'yes')
    cfg           = [];
    cfg.bpfilter  = proc_vars.eog_bp_yn;
    cfg.bpfreq    = proc_vars.eog_bp_freq;
    cfg.bpfiltord = proc_vars.eog_bp_filtord;
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

% Identify most correlated components
avg_eog_ic_corr = mean(eog_ic_corr,3);
heog_ics = find(abs(avg_eog_ic_corr(1,:))>proc_vars.eog_ic_corr_cut);
veog_ics = find(abs(avg_eog_ic_corr(2,:))>proc_vars.eog_ic_corr_cut);
if any([isempty(heog_ics), isempty(veog_ics)])
    error('No EOG ICs found!');
end

%% Generate Figures
if gen_figs 
    % Plot EOG-ICA Correlations
    figure('Visible',fig_vis); hold on;
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
    
    % ICA Power Spectrum, Variance, Topo
    cfg = [];
    cfg.layout   = 'biosemi64.lay';
    cfg.channel  = 'all';
    cfg.path     = SBJ_vars.dirs.proc;
    cfg.prefix   = 'ICA';
    cfg.viewmode = 'component';
    cfg.fig_vis  = fig_vis;
    fn_icabrowser_modified(SBJ, cfg, ica);
    
    % Plot IC single trial stacks + ERPs
    stacked_data_icacomp(SBJ, SBJ_vars.dirs.raw_filename, SBJ_vars.dirs.proc, proc_id, ica);
    
    % Plot IC in ft_databrowser
    if fig_vis
        cfg = [];
        cfg.layout   = 'biosemi64.lay';
        cfg.channel  = 'all';
        cfg.viewmode = 'component';
        ft_databrowser(cfg, ica);
    end
end

%% Save Data
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_clean02a_' proc_id '.mat'];
save(clean_data_fname, '-v7.3', 'trials','cfg_trl','ica', 'heog_ics', 'veog_ics','eog_trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav02a_' proc_id '_clean.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv','bhv_fields');

end
