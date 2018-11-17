%% Load the data
%data_dir = '/Users/SCS22/Desktop/Pilot2_ICATesting/pilot01/'
%data_filename = [data_dir 'ICAsorted_pilot1_take2'];
%data_out_filename = [data_dir 'cleaned_pilot01'];
%load(data_filename, 'eeg', 'eog_bp', 'icaunmixing', 'icatopolabel');

function data_clean(SBJ)
addpath(genpath('/Users/SCS22/Desktop/Knight_Lab/Preprocessed_Work'));
data_in_filename = [SBJ, 'ica+preprocessing']
data_out_filename = [SBJ, 'cleaned']
load(data_in_filename, 'eeg', 'eog_bp', 'icaunmixing', 'icatopolabel');
%% REJECTVISUAL
% Summary

cfg = [];
cfg.method   = 'summary';  % 'summary' for trials+channels; 'channel' for individual trials
clean   = ft_rejectvisual(cfg, eeg);
rejected(1) = {input('specify componenets for rejection in an array: eeg summary method')}
cfg = [];
cfg.method   = 'summary';   % 'summary' for trials+channels; 'channel' for individual trials
clean   = ft_rejectvisual(cfg, eog_bp);
rejected(2) = {input('specify componenets for rejection in an array: eog summary method')}
% Trials
cfg = [];
cfg.method   = 'channel';   % 'summary' for trials+channels; 'channel' for individual trials
clean   = ft_rejectvisual(cfg, eeg);
rejected(3) = {input('specify componenets for rejection in an array: eeg channel method')}
cfg = [];
cfg.method   = 'channel';   % 'summary' for trials+channels; 'channel' for individual trials
clean   = ft_rejectvisual(cfg, eog_bp);
rejected(4) = {input('specify componenets for rejection in an array: eog channel method')}
%% Rejected
eegsummary = rejected{1}
eogsummary = rejected{2}
eegchannel = rejected{3}
eegchannel = rejected{4}
bad_trials = input('trials to remove later in an array format')
bad_channels = input('channels to remove later in a cell array format')

save([data_out_filename 'rejectedtrials'] , 'rejected', 'bad_trials', 'bad_channels');
% Summary vs. Trial Comparison
% compare the lists of trials rejected based on eeg  summary vs. trial, eog
% summary vs. trial, and then all eeg vs. all eog
%% Examine the components
% Rebuild the components
cfg           = [];
cfg.unmixing  = icaunmixing;
cfg.topolabel = icatopolabel;
ica           = ft_componentanalysis(cfg, eeg);

%% eog vs ica Correlation
% low pass filter at 8 hz right before correlatin (dont do before ica) for
% eog_bp and ica
eog_corr = cat(2, [1:numel(ica.topolabel)]', zeros([length(ica.topolabel), 2]));
for x = 1:length(ica.trial) %number trials
    for y = 1:numel(ica.label)
        for z = 2: 3
        % CWH: why do we start at 2, AKA what's in the first column? SS: the
        % labels of component numbers
            temp = corrcoef(eog_bp.trial{1,x}((z-1),:), ica.trial{1,x}(y,:)); %components
            eog_corr(y, z) = temp(1,2);
    end
    end
end

%% Select Top components
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

save(data_out_filename, 'eog_corr');



%% Artifacts1
figure
cfg = [];
cfg.component = top_comps(:);       % specify the component(s) that should be plotted
cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, ica)
cfg.channel = top_comps(:); 
cfg.viewmode  = 'component';
ft_databrowser(cfg, ica);
%% Artifacts
figure
cfg = [];
cfg.component = 1:32;       % specify the component(s) that should be plotted
cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, ica)
figure
cfg = [];
cfg.component = 33:numel(ica.topolabel);       % specify the component(s) that should be plotted
cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, ica)
savefig(data_out_filename)

%% View the components with time seires
cfg.viewmode  = 'component';
cfg.layout    = 'biosemi64.lay'; % specify the layout file that should be used for plotting
cfg.fontsize  = 1/numel(ica.topolabel);
cfg.linewidth = 1;
ft_databrowser(cfg, ica);
%savefig(data_out_filename)

%% ICA Component rejection
cfg = []
rejcomp = input('specify components for rejection in an array format');
cfg.component = rejcomp;
data = ft_rejectcomponent(cfg,ica, data);
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
cfg.channel = {'all', '-AF8'}; %% here how to set the difference between all the channels and bad channels
data=ft_selectdata(cfg, data);
% based on the ft_rejectvisual, toss the bad trials from your
% ICA-reconstructed dataset (which I think should still have all trials...)

%% Save the final clean data
save(data_out_filename, 'data')
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





