%% preprocessing --ERP practice
%data_dir = '/Users/SCS22/Desktop/Pilot2_ICATesting/pilot01/'
%data_filename = [data_dir 'Pilot01.bdf'];
%data_out_filename = [data_dir 'ICAsorted_pilot1_take2'];
%top_comp_cut = 0.1;
%replacement_channels = {{'1-CP2', '1-P2', '1-P8'}, {'1-EX2', '1-EX3', '1-EX6'}};
%ears = {'1-EXG1', '1-EXG2'};
%horiz = {'1-EXG3', '1-EXG4'};
%vert = {'1-EXG5', '1-Fp2'};
%addpath('/Users/SCS22/Documents/MATLAB/fieldtrip/');
%ft_defaults
function SBJ01_import_data(SBJ)

addpath(genpath('/Users/SCS22/Desktop/Preprocessing_Work/'));
addpath('/Users/SCS22/Documents/MATLAB/fieldtrip/');
ft_defaults

%% Load and preprocess the data
SBJ_vars_cmd = ['run ' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
data_out_filename = [SBJ, 'ica_preprocessing']
%% Prepare to Import Data
% eval(strcat(scripts_dir,'SBJName_vars.m'));

% Define trials to prepare to cut into trials
cfg=[];
cfg.dataset             = SBJ_vars.dirs.raw_filename;
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = 2;
cfg.trialdef.prestim    = -.3;
cfg.trialdef.poststim   = 1;
cfg.trialfun            = 'tt_trialfun';
cfg = ft_definetrial(cfg);

%% referencing
cfg.reref = 'yes';
cfg.implicitref   = []; %REF is implicit, what is 0d to, RM is the other one, do average reference with the mastoid
cfg.refchannel    = SBJ_vars.ears;   % L and R ear lobe (didn't record mastoids)

% Baseline-correction options- Not doing right now
%cfg.demean          = 'yes'; %baseline correct data before you resample?
%cfg.baselinewindow  = [-0.25 -0.05]; %set baseline to -.2 to 0 seconds, then correct it before you resample

% Filtering
cfg.lpfilter = 'no';
cfg.hpfilter = 'no';
cfg.bpfilter = 'yes';
%low pass filter at 8 hz
cfg.bpfreq   = [0.1 200]; 
data = ft_preprocessing(cfg);
%% Cut to analysis_time
if ~isempty(SBJ_vars.analysis_time)
    for epoch_ix = 1:length(SBJ_vars.analysis_time)
        epoch_len{epoch_ix} = diff(SBJ_vars.analysis_time{epoch_ix});
        cfg_trim = [];
        cfg_trim.latency = SBJ_vars.analysis_time{epoch_ix};
        
        data_pieces{epoch_ix} = ft_selectdata(cfg_trim, data);
        evnt_pieces{epoch_ix} = ft_selectdata(cfg_trim, evnt);
        if ~isempty(SBJ_vars.ch_lab.eeg)
            eeg_pieces{epoch_ix}  = ft_selectdata(cfg_trim, eeg);
        end
        if ~isempty(SBJ_vars.ch_lab.eog)
            eog_pieces{epoch_ix}  = ft_selectdata(cfg_trim, eog);
        end
    end
    % Stitch them back together
    data = data_pieces{1};
    data.time{1} = data.time{1}-SBJ_vars.analysis_time{1}(1);
    evnt = evnt_pieces{1};
    evnt.time{1} = evnt.time{1}-SBJ_vars.analysis_time{1}(1);
    if ~isempty(SBJ_vars.ch_lab.eeg)
        eeg = eeg_pieces{1};
        eeg.time{1} = eeg.time{1}-SBJ_vars.analysis_time{1}(1);
    end
    if ~isempty(SBJ_vars.ch_lab.eog)
        eog = eog_pieces{1};
        eog.time{1} = eog.time{1}-SBJ_vars.analysis_time{1}(1);
    end
    if length(SBJ_vars.analysis_time)>1
        for epoch_ix = 2:length(SBJ_vars.analysis_time)
            data.trial{1} = horzcat(data.trial{1},data_pieces{epoch_ix}.trial{1});
            data.time{1} = horzcat(data.time{1},data_pieces{epoch_ix}.time{1}-...
                SBJ_vars.analysis_time{epoch_ix}(1)+data.time{1}(end)+data.time{1}(2));
            
            evnt.trial{1} = horzcat(evnt.trial{1},evnt_pieces{epoch_ix}.trial{1});
            evnt.time{1} = horzcat(evnt.time{1},evnt_pieces{epoch_ix}.time{1}-...
                SBJ_vars.analysis_time{epoch_ix}(1)+evnt.time{1}(end)+evnt.time{1}(2));
            
            if ~isempty(SBJ_vars.ch_lab.eeg)
                eeg.trial{1} = horzcat(eeg.trial{1},eeg_pieces{epoch_ix}.trial{1});
                eeg.time{1} = horzcat(eeg.time{1},eeg_pieces{epoch_ix}.time{1}-...
                    SBJ_vars.analysis_time{epoch_ix}(1)+eeg.time{1}(end)+eeg.time{1}(2));
            end
            if ~isempty(SBJ_vars.ch_lab.eog)
                eog.trial{1} = horzcat(eog.trial{1},eog_pieces{epoch_ix}.trial{1});
                eog.time{1} = horzcat(eog.time{1},eog_pieces{epoch_ix}.time{1}-...
                    SBJ_vars.analysis_time{epoch_ix}(1)+eog.time{1}(end)+eog.time{1}(2));
            end
        end
    end
end

%% Resampling
cfg=[];
cfg.resamplefs = 250;
data = ft_resampledata(cfg, data);

cfg = [];
badindex = match_str(data.label, SBJ_vars.replacement_channels{1});
cfg.channel = setdiff([1:numel(data.label)], badindex);
%cfg.channel = setdiff(data.label,replacementchannels{1});
data = ft_selectdata(cfg, data);
%try to remove the replacement channels here by changing channel, wasnt
%workign with strings
%% Get rid of 1s prefix
for i = 1: numel(data.label)
   data.label{i} = data.label{i}(3:end);
   for x = 1: numel(SBJ_vars.replacement_channels{1})
        if (strcmp(SBJ_vars.replacement_channels{2}{x},data.label{i}))
           data.label{i} = SBJ_vars.replacement_channels{1}{x}; % replaces the label of the externals with the channel they represent
        end
   end
end
%% Extract EOGs and Bipolar Rereference (based off Arjen's)
% Extract EOG channels
cfg = [];
cfg.channel = ft_channelselection({'EXG*', 'Fp2'}, data.label);  % before was using ft_channelselection
eog = ft_selectdata(cfg, data);

%% Prepare the bipolar version
cfg.channel = SBJ_vars.ears; %%had some trouble with this line working too
eog_bp = eog;
horizindex1 = match_str(eog.label, SBJ_vars.horiz{1});
vertindex1 = match_str(eog.label, SBJ_vars.vert{1});
horizindex2 = match_str(eog.label, SBJ_vars.horiz{2});
vertindex2 = match_str(eog.label, SBJ_vars.vert{2});
for t = 1:length(eog.trial)
    %!!! get rid fo the ear channels!
    eog_bp.trial{1, t}(1, :) = eog.trial{1,t}(horizindex1, :) - eog.trial{1,t}(horizindex2, :); %% again, how to select strings instead of indexing into cell array?? to pick out horiz(1) - horiz(2)
    eog_bp.trial{1, t}(2, :) = eog.trial{1,t}(vertindex1, :) - eog.trial{1,t}(vertindex2, :); %vertical
end

%% Select EEG Data
cfg = [];
cfg.channel = {'all', '-EXG*', '-atus'}; % leave in FP2
eeg = ft_selectdata(cfg, data);

%% ICA
cfg        = [];
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
icomp = ft_componentanalysis(cfg, eeg);

%% save the ICA unmixing and topolabel, eeg, and eog_bp
icaunmixing = icomp.unmixing;
icatopolabel = icomp.topolabel;
save(data_out_filename, 'icaunmixing', 'icatopolabel', 'eeg', 'eog_bp');


