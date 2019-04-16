function EEG02a_artifact_rejection(SBJ, proc_id, plt_id,genfigs, figon)

if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);
plot_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plot_vars_cmd);
data_fname = [SBJ_vars.dirs.preproc SBJ '_preproc_' proc_id '.mat'];
load(data_fname);

%% Import behavioral data
%   Total_Trial,Block,Condition,Hit,RT,Timestamp,Tolerance,Trial,Score,ITI,ITI type
fprintf('\tReading behavioral csv file\n');
bhv_file = fopen([SBJ_vars.dirs.events SBJ '_behav.csv'], 'r');
bhv_fields = textscan(bhv_file, '%s %s %s %s %s %s %s %s %s %s %s', 1, 'Delimiter', ',');
bhv_data = textscan(bhv_file, '%d %d %s %d %f %f %f %d %d %f %f',...
    'Delimiter', ',', 'HeaderLines', 1);
fclose(bhv_file);

for f_ix = 1:numel(bhv_fields)
    bhv_fields{f_ix} = strrep(bhv_fields{f_ix}{1},' ','_');
    bhv.(bhv_fields{f_ix}) = bhv_data{f_ix};
end

%% Cut into trials
% Must segment before downsampling because trigger channel read from
% original file
%Need to recut trials on updated data with the nans
cfg = [];
cfg.dataset             = SBJ_vars.dirs.raw_filename;
cfg.trialdef.eventtype  = 'STATUS';%SBJ_vars.ch_lab.trigger;
cfg.trialdef.eventvalue = proc_vars.event_code;        % feedback cocde
cfg.trialdef.prestim    = proc_vars.trial_lim_s(1);
cfg.trialdef.poststim   = proc_vars.trial_lim_s(2);
cfg.trialfun            = 'tt_trialfun';%'ft_trialfun_general';%
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

trials = ft_redefinetrial(cfg_trl,data);
eog_trials = ft_redefinetrial(cfg_trl,eog);
%% Exclude bad_trials
load([SBJ_vars.dirs.events SBJ '_raw_bad_epochs.mat']);

%Find epochs that overlap with bad_epochs from raw visual inspection
if ~isempty(bad_epochs)
    bad_raw_trials = fn_find_trials_overlap_epochs(bad_epochs,1:size(data.trial{1},2),...
        event_onsets,proc_vars.trial_lim_s*data.fsample);    
else
    bad_raw_trials = [];
end

% Exclude bad trials
cfg_trl.trl(bad_raw_trials,:) = [];
trials = ft_redefinetrial(cfg_trl,trials);
eog_trials = ft_redefinetrial(cfg_trl,eog);


for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(bad_raw_trials) = [];
end

% Identify and exclude training and bad raw visual trials
training_trial_ix = find(bhv.Block==-1); % returns index of trial number -- one indexed
response_time_low = find(bhv.RT <=proc_vars.rt_bounds(1)); % returns index of trial number -- one indexed
response_time_high = find(bhv.RT>= proc_vars.rt_bounds(2)); % returns index of trial number - one indexed
exclude_trials = unique(vertcat(training_trial_ix,response_time_low, response_time_high));
bad_trials = exclude_trials';

for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(exclude_trials) = [];
end

% Resegment trials
cfg_trl.trl(bad_trials,:) = [];
trials = ft_redefinetrial(cfg_trl,trials);
eog_trials = ft_redefinetrial(cfg_trl,eog);



%% EOG vs. ICA Correlation
% Rebuild the components
cfg           = [];
cfg.unmixing  = icaunmixing;
cfg.topolabel = icatopolabel;
ica           = ft_componentanalysis(cfg, trials);

% Filter EOG
cfg           = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq    = [1 15];  % from ft_rejectvisual help page
cfg.bpfiltord = 4;       % from ft_rejectvisual help page
eog_bp = ft_preprocessing(cfg,eog_trials);

eog_ic_corr = zeros([numel(eog_bp.label), numel(ica.topolabel), numel(eog_bp.trial)]);
    for eog_ix = 1:2
        for ic_ix = 1:numel(ica.label)
            for t_ix = 1:length(ica.trial) %number trials
                temp = corrcoef(eog_bp.trial{t_ix}(eog_ix,:), ica.trial{t_ix}(ic_ix,:));
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
%% Downsample
%Just looked over script a final time! This is one thing I wanted to check
%that its in the right spot that I forgot to mention -- I wanted it to go
%after filtering and trial cutting, but wasn't sure if this was too late!
%(April 13,2019)
if strcmp(proc_vars.resample_yn,'yes')
      cfg=[];
    cfg.resamplefs = proc_vars.resample_freq;
    trials = ft_resampledata(cfg, trials);
    eog_trials = ft_resampledata(cfg, eog_trials);
end

if genfigs 
    %% Graph EOG Correlations
    figure;
    subplot(1,2,1);
    histogram(avg_eog_ic_corr(1,:),50);
    title('HEOG ICA Correlations');
    xlabel('Correlation (r)');
    subplot(1,2,2);
    histogram(avg_eog_ic_corr(2,:),50);
    title('VEOG ICA Correlations');
    xlabel('Correlation (r)');
    eog_ic_fname = [SBJ_vars.dirs.preproc SBJ '_eog_ic_corr_' proc_id '.png'];
    saveas(gcf,eog_ic_fname);

    fprintf('==================== EOG vs. ICA Correlations =========================\n');
    fprintf(['%i ICs correlate with HEOG, r = ' repmat('%.02f ',[1 numel(heog_ics)]) '\n'],...
            numel(heog_ics),avg_eog_ic_corr(1,heog_ics));
    fprintf(['%i ICs correlate with VEOG, r = ' repmat('%.02f ',[1 numel(veog_ics)]) '\n'],...
            numel(veog_ics),avg_eog_ic_corr(2,veog_ics));
    fprintf('=======================================================================\n');
    %% POWER SPECTRUM, VARIANCE, TOPO
    cfg = [];
    cfg.layout = 'biosemi64.lay';
    cfg.channel = 'all'
    cfg.path = SBJ_vars.dirs.proc
    cfg.prefix = 'ICA';
    cfg.viewmode = 'component';
    ft_icabrowser_modified(SBJ, cfg, ica);
    %% stacked data
    stacked_data_icacomp(SBJ, SBJ_vars.dirs.raw_filename, SBJ_vars.dirs.proc, proc_id, ica);
end

    clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_clean02a_' proc_id '.mat'];
    save(clean_data_fname, '-v7.3', 'trials','cfg_trl','ica', 'heog_ics', 'veog_ics','eog_trials');

    clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav02a_' proc_id '_clean.mat'];
    save(clean_bhv_fname, '-v7.3', 'bhv','bhv_fields');
    
if figon
    load([root_dir 'PRJ_Error_eeg/scripts/utils/cfg_plot_eeg.mat']);
    cfg = [];
    cfg.layout    = 'biosemi64.lay';
    cfg.channel = 'all'; 
    cfg.viewmode  = 'component';
    ft_databrowser(cfg, ica);
end
