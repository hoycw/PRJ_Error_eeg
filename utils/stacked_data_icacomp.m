
function stacked_data_icacomp(SBJ, SBJ_dir, SBJ_dirproc, proc_id, data)
%% Data Preparation
% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';ft_dir='/Users/SCS22/Documents/MATLAB/fieldtrip/';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load the data

%% Import behavioral data
%   Total_Trial,Block,Condition,Hit,RT,Timestamp,Tolerance,Trial,Score,ITI,ITI type

events = ft_read_event(SBJ_dir);
eventvalues = extractfield(events, 'value');
eventsamples = extractfield(events, 'sample');
if numel(eventvalues) ~= numel(eventsamples)
   difference = numel(eventsamples)-numel(eventvalues);
   eventsamples = eventsamples(difference+1:end);
end
stimulus_ix = find(eventvalues == 1);
feedback_ix = find(eventvalues == 2);
stimulus_onset = eventsamples(stimulus_ix);
feedback_onset = eventsamples(feedback_ix);
stim_onset = find(data.time{1,1} == 0);
stim_onset_array(1:numel(stimulus_onset)) = stim_onset;
[~,fb_onset] = min(abs(1.8-data.time{1,1}));
fb_onset_array(1:numel(feedback_onset)) = fb_onset;
[~,one_sec] = min(abs(1 - data.time{1,1}));
one_sec_array(1:numel(feedback_onset)) = one_sec;
%for x = 1: numel(feedback_onset);
    %trial_sample = feedback_onset(x)/numel(data.time{1,x}(1,:))

cfg = [];
cfg.channel = 'all';
erps = ft_timelockanalysis(cfg, data);
            % keep manual screen position - better in dual monitor settings

        for comp_ix = 1:numel(data.label)
            fig_name = [SBJ data.label{comp_ix}];
            f=figure('units','normalized','Name', fig_name, 'Visible', 'off');
            subplot('Position', [0.1, 0.35, 0.8, 0.6]); 
            clims = NaN([1 2]);
            comp_data = NaN([numel(data.trial) numel(data.time{1})]);
            for t_ix = 1:numel(data.trial)
                comp_data(t_ix, :) = squeeze(data.trial{t_ix}(comp_ix,:));
            end
            clims(1,1) = prctile(comp_data(:),5);
            clims(1,2) = prctile(comp_data(:),95);
            clims = [min(clims(:,1)) max(clims(:,2))];
            imagesc(comp_data);
            set(gca,'YDir','normal');
            ax = gca;
            ax.XTick = [1, 257, 513, 769, 1025, 1281, 1537, 1793, 2049, 2305, 2561, 2817, 3073];
            ax.XTickLabel = [-.25, 0, 0.25, .5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75];
            y = linspace(1,numel(stimulus_onset), numel(stimulus_onset));
            y2 = linspace(1,numel(feedback_onset), numel(feedback_onset));
            stim_line = line(stim_onset_array,y,...
                    'LineWidth',2,'Color','k');
            feedback_line = line(fb_onset_array,y2,...
                    'LineWidth',2,'Color','k');
            one_line = line(one_sec_array, y2, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--');
            colorbar
            subplot('Position', [0.1, 0.1, 0.70, 0.15]); 
            plot(erps.time,erps.avg(comp_ix,:),'k');
            xlim([-.25 2.8])
            comp_stack_name = [SBJ_dirproc SBJ data.label{comp_ix,1} '_fig2' '.png'];
            saveas(gcf,comp_stack_name);      
        end
