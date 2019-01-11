function EEG02_artifact_rejection(SBJ, proc_id, dorejectvisual)

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath(genpath([root_dir 'PRJ_Error_eeg/scripts/']));
addpath(ft_dir);
ft_defaults

%% Load the data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_proc_vars.m'];
eval(proc_vars_cmd);

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

%% Exclude training and bad visual trials
% Check that behavioral and EEG event triggers line up
if numel(bhv.Trial)~=numel(event_onsets)
    error(['Mismatch in behavioral and neural trial counts: ' num2str(numel(bhv.Trial))...
        ' behavioral; ' num2str(numel(event_onsets)) ' neural']);
end

% Identify and exclude training and bad raw visual trials
training_trial_n = find(bhv.Block==-1);
exclude_trials = unique(vertcat(training_trial_n,bad_raw_trials));
for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(exclude_trials) = [];
end

% Resegment trials
cfg_trl.trl(exclude_trials,:) = [];
trials = ft_redefinetrial(cfg_trl,data);
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

% Correlation between both EOGs and all ICs for all trials
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

% IC rejection
cfg = [];
cfg.component = [heog_ics veog_ics];
clean_trials = ft_rejectcomponent(cfg, ica, trials);

% Print the results
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

% View ICs with EOG correlations
cfg = [];
cfg.layout    = 'biosemi64.lay';
cfg.channel = [heog_ics veog_ics]; 
cfg.viewmode  = 'component';
ft_databrowser(cfg, ica);

%% Visual Trial Rejection
if dorejectvisual
    cfg = [];
    cfg.method = 'summary';  % 'summary' for trials+channels; 'channel' for individual trials
    clean_summ = ft_rejectvisual(cfg, clean_trials);
    
    % Trials
    cfg = [];
    cfg.method  = 'channel';   % 'summary' for trials+channels; 'channel' for individual trials
    clean_trial = ft_rejectvisual(cfg, clean_trials);
    
    % Report channels and trials identified above in SBJ_vars, then re-run
    exit;
end

%% Clean up and save data
% Get bad trials and channels from SBJ_vars
%   NOTE: these are indices into the post-raw rejection trial list
bad_trials = SBJ_vars.trial_reject_ix;
bad_ch_neg = fn_ch_lab_negate(SBJ_vars.ch_lab.bad);

% Remove from ICA cleaned data
cfgs = [];
cfgs.trials = setdiff(1:numel(clean_trials.trial),bad_trials);
cfgs.channel = {'all',bad_ch_neg{:}};
trials = ft_selectdata(cfgs,clean_trials);

% Remove from behavioral
for f_ix = 1:numel(bhv_fields)
    bhv.(bhv_fields{f_ix})(bad_trials) = [];
end
    
% Save outputs
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_clean_' proc_id '.mat'];
save(clean_data_fname, '-v7.3', 'trials');

clean_bhv_fname = [SBJ_vars.dirs.events SBJ '_behav_' proc_id '_clean.mat'];
save(clean_bhv_fname, '-v7.3', 'bhv');

