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
data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '.mat'];
load(data_fname);

%% Import behavioral data
%   Total_Trial,Block,Condition,Hit,RT,Timestamp,Tolerance,Trial,Score,ITI,ITI type
fprintf('\tReading behavioral csv file\n');
bhv_file = fopen([SBJ_vars.dirs.events SBJ '_behav.csv'], 'r');
bhv_fields = textscan(bhv_file, '%s %s %s %s %s %s %s %s %s %s %s', 1, 'Delimiter', ',');
bhv_data = textscan(bhv_file, '%d %d %s %d %f %f %f %d %d %f %f',...
    'Delimiter', ',', 'HeaderLines', 1);
fclose(bhv_file);

for t_ix = 1:numel(bhv_data{1})
    for f_ix = 1:numel(bhv_fields)
        bhv.(strrep(bhv_fields{f_ix}{1},' ','_')) = bhv_data{f_ix};
    end
end

%% REJECT VISUAL
% Match BHV vs EEG trials
% Exclude example and training data
% Exclude outlier RTs?
if (isempty(SBJ_vars.trial_reject_ix) || dorejectvisual)
    cfg = [];
    cfg.method   = 'summary';  % 'summary' for trials+channels; 'channel' for individual trials
    clean   = ft_rejectvisual(cfg, trials);
    % specify trials and/or channels for rejection in an array:
    
    cfg = [];
    cfg.method   = 'summary';   % 'summary' for trials+channels; 'channel' for individual trials
    clean   = ft_rejectvisual(cfg, eog);
    rejected(2) = {input('specify componenets for rejection in an array: eog summary method')}
    % Trials
    cfg = [];
    cfg.method   = 'channel';   % 'summary' for trials+channels; 'channel' for individual trials
    clean   = ft_rejectvisual(cfg, trials);
    rejected(3) = {input('specify componenets for rejection in an array: eeg channel method')}
    cfg = [];
    cfg.method   = 'channel';   % 'summary' for trials+channels; 'channel' for individual trials
    clean   = ft_rejectvisual(cfg, eog);
    rejected(4) = {input('specify componenets for rejection in an array: eog channel method')}
    eegsummary = rejected{1}
    eogsummary = rejected{2}
    eegchannel = rejected{3}
    eegchannel = rejected{4}
    bad_trials = input('trials to remove later in an array format')
    bad_channels = input('channels to remove later in a cell array format')
else
    bad_trials = SBJ_vars.trial_reject_ix;
    bad_channels = SBJ_vars.channels_reject_ix;
end
!!! saving out the rejected trials + channels, etc.
save([clean_fname 'rejectedcomponents'] , 'bad_trials', 'bad_channels');
% Summary vs. Trial Comparison
% compare the lists of trials rejected based on eeg  summary vs. trial, eog
% summary vs. trial, and then all eeg vs. all eog
%% Examine the components
% Rebuild the components
cfg           = [];
cfg.unmixing  = icaunmixing;
cfg.topolabel = icatopolabel;
ica           = ft_componentanalysis(cfg, trials);

%% eog vs ica Correlation
% low pass filter at 8 hz right before correlatin (dont do before ica) for
% eog_bp and ica
cfg           = [];
cfg.bpfilter  = 'yes';
cfg.bpfreq    = [1 15];  % from ft_rejectvisual help page
cfg.bpfiltord = 4;       % from ft_rejectvisual help page
eog_bp = ft_preprocessing(cfg,eog_trials);
eog_corr = cat(2, [1:numel(ica.topolabel)]', zeros([length(ica.topolabel), 2]));
for x = 1:length(ica.trial) %number trials
    for y = 1:numel(ica.label)
        for z = 2:3 % 
        % CWH: why do we start at 2, AKA what's in the first column? SS: the
        % labels of component numbers
            temp = corrcoef(eog_trials.trial{1,x}((z-1),:), ica.trial{1,x}(y,:)); %components
            eog_corr(y, z) = temp(1,2);
        end
    end
end

%% Select Top components based on correlation values
!!! get right top_cut from proc_vars
num_top_comp = ceil(size(eog_corr,1)*SBJ_vars.top_comp_cut);
top_comps = zeros(num_top_comp,6);
[correlation, top_horiz_ix] = sort(abs(eog_corr(:,2)),'descend');
top_comps(:,1) = eog_corr(top_horiz_ix(1:num_top_comp),1);
top_comps(:,2) = eog_corr(top_horiz_ix(1:num_top_comp),2);
top_comps(:,3) = eog_corr(top_horiz_ix(1:num_top_comp),3);
[correlation, top_vert_ix] = sort(abs(eog_corr(:,3)),'descend');
top_comps(:,4) = eog_corr(top_vert_ix(1:num_top_comp),1);
top_comps(:,5) = eog_corr(top_vert_ix(1:num_top_comp),2);
top_comps(:,6) = eog_corr(top_vert_ix(1:num_top_comp),3);

top_comps
% Print out some useful number on the results (top compnoetn #'s, their
% correlation values, etc.)
% this prints out the correlation with the horizontal in column 2/5 and the
% vertical in column 3/6 -- not sure if there is a way to label this
!!! clean_fname = clean' datestr(now,'mm-dd-yyyy HH-MM')

save(clean_fname, 'eog_corr');

%% Viewing ICA components based on top correlations
figure
cfg = [];
cfg.component = top_comps(:, [1 4]);       % specify the component(s) that should be plotted
cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, ica)
!!! get rid of this top (topo-only) part because the ft_databrowser component view also shows the topoplots
cfg.channel = top_comps(:, [1 4]); 
cfg.viewmode  = 'component';
ft_databrowser(cfg, ica);
% %% Artifacts
% figure
% cfg = [];
% cfg.component = 1:32;       % specify the component(s) that should be plotted
% cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';
% ft_topoplotIC(cfg, ica)
% figure
% cfg = [];
% cfg.component = 33:numel(ica.topolabel);       % specify the component(s) that should be plotted
% cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';
% ft_topoplotIC(cfg, ica)
% savefig(data_out_filename_cleaned)

%% View all ICA components with time seires
cfg.viewmode  = 'component';
cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
cfg.fontsize  = 1/numel(ica.topolabel);
cfg.linewidth = 1;
ft_databrowser(cfg, ica);
%savefig(data_out_filename_cleaned)

%% ICA Component rejection
cfg = []
!!! don't use input, jsut manually enter in SBJ_vars and then load that.
rejcomp = input('specify components for rejection in an array format');
cfg.component = rejcomp;
data = ft_rejectcomponent(cfg,ica, trials);
% after selecting a couple components to reject, rebuild the data without them
% If you aren't sure which ones to reject (likely), just pick your top
% couple and do it anyways so you have the code to do it once we meet and
% decide together...

%% Trial rejection
%cfg = [];
%cfg.channel =setdiff('all', bad_channels);
%data3.channel = data2.label(cfg.channel);
cfg = [];
cfg.trial = setdiff([1: numel(data.trial)], bad_trials);
data.trial = data.trial(cfg.trial);
data.time   = data.time(cfg.trial);
cfg=[];
badindex = match_str(data.label, bad_channels);
cfg.channel = setdiff([1:numel(data.label)], badindex); %%here how to set the difference between all the channels and bad channels
data=ft_selectdata(cfg, data);
% based on the ft_rejectvisual, toss the bad trials from your
% ICA-reconstructed dataset (which I think should still have all trials...)
cfg.viewmode = 'vertical';
browsed_data_raw = ft_databrowser(cfg, data);
!!!ica_fname = 'Pilot02ICASorted' datestr(now,'mm-dd-yyyy HH-MM')
save(data_out_filename_ICA, 'browsed_data_raw');
artifacts_data_raw = browsed_data_raw.artfctdef.visual.artifact;
rawdata = cfg.dataset;
cleandata = ft_rejectartifact(cfg,rawdata);


%% Save the final clean data
save(data_out_filename_cleaned, 'data')
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





