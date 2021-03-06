function SBJ06a_OB_ERP_save_p2p(SBJ_id,proc_id,feat_id,varargin)
%% Save amplitude and latency of peak-to-peak N2s in oddball (OB) task
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load and average single-trial ERP within condition
%   Find positive and negative peaks:
%       Must find positive and negative peaks in their respective windows
%       Positive peak must precede negative peak
%   Compute amplitude differences:
%       For multiple peak pairs, selects largest amplitude difference
%       Rejects pairs with amplitude difference in wrong direction (negative > positive)
%   Plot peak detection errors
%       Optional: Plot ERPs and detected peaks
%   Save amplitudes and latencies
% INPUTS:
%   SBJ_id [str]  - ID of subject list for group
%   proc_id [str] - ID of oddball preprocessing pipeline
%   feat_id [str] - ID of the feature extraction parameters to use
%       ft.an_id  = 'ERP_all_S2t1_dm2t0_fl05t20'
%       ft.grp_id = conditions to extract features (likely 'Tar','Odd','rare')
%       ft.measure = 'p2p', difference in amplitude of preceding + and following - peaks
%   varargin:
%       plot_feat [0/1] - binary flag for plotting each ERP with peaks overlay
%           default: 1
%       plot_peaks [0/1] - binary flag for peak detection errors and summary plot of both peaks
%           default: 1
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)
%   ft_amp [float] - amplitudes of [positive, negative] features
%   ft_times [float] - times (in sec) of [positive, negative] features
%   pk_amp [float] - amplitudes of both peaks per SBJ
%   pk_times [float] - latency (in sec) of both peaks per SBJ
%   miss_fts [boolean] - array indicating if features couldn't be extracted per SBJ
%   bad_pks [boolean] - array indicating if peaks didn't fit P2P per SBJ

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Handle Variable Inputs & Defaults
if ~isempty(varargin)
    for v = 1:2:numel(varargin)
        if strcmp(varargin{v},'fig_vis') && ischar(varargin{v+1})
            fig_vis = varargin{v+1};
        elseif strcmp(varargin{v},'fig_ftype') && ischar(varargin{v+1})
            fig_ftype = varargin{v+1};
        elseif strcmp(varargin{v},'plot_peaks')
            plot_feat = varargin{v+1};
        elseif strcmp(varargin{v},'plot_feat')
            plot_feat = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');    fig_vis = 'on'; end
if ~exist('fig_ftype','var');  fig_ftype = 'png'; end
if ~exist('save_fig','var');   save_fig = 1; end
if ~exist('plot_feat','var');  plot_feat = 1; end
if ~exist('plot_peaks','var'); plot_peaks = 1; end

%% Load Data 
if ~contains(proc_id,'odd'); error('proc_id must be for oddball task!'); end
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' feat_id '_vars.m'];
eval(feat_vars_cmd);
if ~strcmp(ft.measure,'p2p'); error('This script is only for peak-to-peak!'); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, ~, cond_colors, ~, ~] = fn_condition_label_styles(ft.grp_id);

fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' ft.an_id '/OB_feat_ts/'];
if ~exist(fig_dir,'dir') && save_fig
    mkdir(fig_dir);
end

%% Load Behavior and Select Conditions
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
for s = 1:numel(SBJs)
    % Load data
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    bhvs{s} = tmp.bhv;
    
    % Select Conditions of Interest
    full_cond_idx{s} = fn_condition_index(cond_lab, bhvs{s});
    bhv_fields = fieldnames(bhvs{s});
    orig_n_trials = numel(bhvs{s}.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhvs{s}.(bhv_fields{f_ix}))==orig_n_trials
            bhvs{s}.(bhv_fields{f_ix}) = bhvs{s}.(bhv_fields{f_ix})(full_cond_idx{s}~=0);
        end
    end
    n_trials(s) = numel(bhvs{s}.trl_n);
    
    clear tmp
end

%% Load Data and Compute ERPs
cfgs  = [];
% cfgs.latency = [min(ft.lim(:)) max(ft.lim(:))];
cfgs.channel = unique(ft.chan);
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',ft.an_id,'.mat'],'roi');
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    ft_roi = ft_selectdata(cfgs, roi);
    
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = ft_roi.time{1};
        ch_list  = ft_roi.label;
        erps     = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
    end
    
    % Compute ERPs within condition
    cond_idx = fn_condition_index(cond_lab, bhvs{s});
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for trl_ix = 1:numel(cond_trial_ix)
            trials(:,trl_ix,:) = ft_roi.trial{cond_trial_ix(trl_ix)};
        end
        erps(cond_ix,s,:,:) = mean(trials,2);
    end
    
    clear tmp roi ft_roi trials cond_idx cond_trial_ix
end

%% Compute Peak Amplitude and Latency
ft_amp   = nan([numel(SBJs) numel(ft.name)]);   % P2P amplitude
% ft_times = nan([numel(SBJs) numel(ft.name)]); % this doesn't make sense...
pk_amp   = nan([numel(SBJs) numel(ft.name) 2]); % individual peak amplitudes
pk_times = nan([numel(SBJs) numel(ft.name) 2]); % individual peak latencies
miss_fts = false([numel(SBJs) numel(ft.name)]);
bad_pks  = false([numel(SBJs) numel(ft.name)]);
for ft_ix = 1:numel(ft.name)
    % Select peak feature parameters
    cond_ix = find(strcmp(cond_lab,ft.cond{ft_ix}));
    ch_ix = find(strcmp(ch_list,ft.chan{ft_ix}));
    
    % Confirm peak signs and order
    if ft.pk_sign{ft_ix}(1)~=1 || ft.pk_sign{ft_ix}(2)~=-1
        error('First peak should be positive, and second should be negative!');
    end
    
    % Find Peak Search Range
    ft_rng = zeros([2 2]);
    ft_time_vec = cell([2 1]);
    for pk_ix = 1:2
        % Match to time vector
        for lim_ix = 1:2
            [~, ft_rng(pk_ix,lim_ix)] = min(abs(time_vec-ft.lim{ft_ix}(pk_ix,lim_ix)));
        end
        ft_time_vec{pk_ix} = time_vec(ft_rng(pk_ix,1):ft_rng(pk_ix,2));
    end
    
    % Peak-to-peak algorithm
    for s = 1:numel(SBJs)
        % Find all possible positive and negative peak amplitudes and latencies
        tmp_amp = cell([2 1]); tmp_times = cell([2 1]);
        for pk_ix = 1:2
            [tmp_amp{pk_ix}, tmp_lat] = findpeaks(ft.pk_sign{ft_ix}(pk_ix)*...
                squeeze(erps(cond_ix,s,ch_ix,ft_rng(pk_ix,1):ft_rng(pk_ix,2))));
            % Convert latencies to time
            if ~isempty(tmp_lat)
                tmp_times{pk_ix} = ft_time_vec{pk_ix}(tmp_lat);
            end
        end
        
        % Check for missing peaks
        if isempty(tmp_amp{1}) || isempty(tmp_amp{2})
            miss_fts(s,ft_ix) = true;
            if isempty(tmp_amp{1})
                sign_str = 'pos';
            else
                sign_str = 'neg';
            end
            fprintf(2,'\tNo %s peak detected for %s %s!\n',...
                sign_str,SBJs{s},cond_lab{cond_ix});
        else
            % Compute peak-to-peak differences using all preceding
            %   positive peaks for each negative peak
            tmp_p2p    = ones(size(tmp_amp{2}))*-999;   % Initialize to tiny value
            p2p_pos_ix = nan(size(tmp_amp{2}));         % Best positive peak per negative peak
            for neg_pk_ix = 1:numel(tmp_amp{2})
                % Only use positive peaks preceding this negative peak
                pre_neg_pos_ix = find(tmp_times{1}<tmp_times{2}(neg_pk_ix));
                for i = 1:numel(pre_neg_pos_ix)
                    % Compare amplitude differential to find greatest
                    %   After sign flip in findpeaks, addition yields subtraction
                    if tmp_amp{1}(pre_neg_pos_ix(i))+tmp_amp{2}(neg_pk_ix) > tmp_p2p(neg_pk_ix)
                        tmp_p2p(neg_pk_ix) = tmp_amp{1}(pre_neg_pos_ix(i))+tmp_amp{2}(neg_pk_ix);
                        p2p_pos_ix(neg_pk_ix) = pre_neg_pos_ix(i);
                    end
                end
            end
            
            % Select pair of peaks with max peak-to-peak difference
            [~, max_p2p_neg_ix] = max(tmp_p2p);
            pk_amp(s,ft_ix,1)   = tmp_amp{1}(p2p_pos_ix(max_p2p_neg_ix));
            pk_times(s,ft_ix,1) = tmp_times{1}(p2p_pos_ix(max_p2p_neg_ix));
            pk_amp(s,ft_ix,2)   = -tmp_amp{2}(max_p2p_neg_ix);  % Flip sign back after findpeaks
            pk_times(s,ft_ix,2) = tmp_times{2}(max_p2p_neg_ix);
        end
        
        % Compute peak-to-peak amplitude difference
        if ~miss_fts(s,ft_ix)
            amp_diff = pk_amp(s,ft_ix,1)-pk_amp(s,ft_ix,2);
            % Toss differences that go in the wrong direction or order
            %   This indicates a complex, multi-peak waveform...
            if amp_diff > 0
                ft_amp(s,ft_ix) = amp_diff;
            else
                bad_pks(s,ft_ix) = true;
                fprintf(2,'\tBad peaks detected for %s %s!\n',...
                    SBJs{s},ft.name{ft_ix});
            end
        end
        
        % Plot Peak Detection Errors (missing or bad peaks)
        if plot_peaks && (any(miss_fts(s,ft_ix)) || bad_pks(s,ft_ix))
            fig_name = ['bad_peaks_' ft.name{ft_ix} '_' SBJs{s} ...
                '_' cond_lab{cond_ix} '_' ch_list{ch_ix}];
            figure('Name',fig_name,'Visible',fig_vis);
            
            % Plot ERP
            plot(time_vec(ft_rng(1,1):ft_rng(2,2)),...
                squeeze(erps(cond_ix,s,ch_ix,ft_rng(1,1):ft_rng(2,2))),...
                'Color',cond_colors{cond_ix});
            hold on;
            
            % Plot Peaks on top
            for pk_ix = 1:2
                [tmp_amp, tmp_lat] = ...
                    findpeaks(ft.pk_sign{ft_ix}(pk_ix)*...
                    squeeze(erps(cond_ix,s,ch_ix,ft_rng(pk_ix,1):ft_rng(pk_ix,2))));
                if ft.pk_sign{ft_ix}(pk_ix)==1; mrkr_color = 'b'; else mrkr_color = 'r'; end
                if ~miss_fts(s,ft_ix)
                    scatter(time_vec(tmp_lat+ft_rng(pk_ix,1)-1),...
                        ft.pk_sign{ft_ix}(pk_ix)*tmp_amp, 'o', mrkr_color);
                end
            end
            
            % Create error message for plot title
            bad_str = [];
            if bad_pks(s,ft_ix)
                bad_str = [bad_str ' bad;'];
            end
            if any(miss_fts(s,ft_ix))
                bad_str = [bad_str ' miss;'];
            else
                % Plot line connecting peaks
                line([pk_times(s,ft_ix,1) pk_times(s,ft_ix,2)],...
                    [pk_amp(s,ft_ix,1) pk_amp(s,ft_ix,2)],...
                    'Color','k');
            end
            
            % Plot parameters
            title([SBJs{s} ' ' cond_lab{cond_ix} ': ' bad_str]);
            xlabel('Time (s)'); ylabel('Amplitude (uV)');
            set(gca,'FontSize',16);
            hold off;
            
            % Save figure
            if save_fig
                bad_fig_dir = [fig_dir 'bad_peak_detection/'];
                if ~exist(bad_fig_dir,'dir'); mkdir(bad_fig_dir); end
                fig_fname = [bad_fig_dir fig_name '.' fig_ftype];
                saveas(gcf,fig_fname);
            else
                pause;
            end
        end
    end
end

%% Plot Features
if plot_feat
    fig_name = [SBJ_id '_' feat_id '_' proc_id];
    if numel(ft.name)<=2
        figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 1 0.5],'Visible',fig_vis);
    else
        figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 1 1],'Visible',fig_vis);
    end
    [n_rowcol,~] = fn_num_subplots(numel(ft.name));
    sbj_colors = distinguishable_colors(numel(SBJs));
    
    for ft_ix = 1:numel(ft.name)
        subplot(n_rowcol(1),n_rowcol(2),ft_ix); hold on;
        cond_ix = strcmp(cond_lab,ft.cond{ft_ix});
        ch_ix   = strcmp(ch_list,ft.chan{ft_ix});
        
        for s = 1:numel(SBJs)
            % Plot ERPs
            plot(time_vec,squeeze(erps(cond_ix,s,ch_ix,:)),'Color',sbj_colors(s,:));
            
            % Plot Peaks
            scatter(pk_times(s,ft_ix,:),pk_amp(s,ft_ix,:), 'o', ...
                'MarkerEdgeColor', sbj_colors(s,:));
            
            % Plot peak-to-peak distance as a line
            line([pk_times(s,ft_ix,1) pk_times(s,ft_ix,2)],...
                [pk_amp(s,ft_ix,1) pk_amp(s,ft_ix,2)],...
                'Color',sbj_colors(s,:),'LineWidth',2);
        end
        
        % Plot parameters
        title([ft.name{ft_ix} ': ' ft.cond{ft_ix} ', ' ft.chan{ft_ix}]);
        xlabel('Time (s)'); ylabel('Amplitude (uV)');
        xlim([ft.plot_lim(1) ft.plot_lim(2)]);
        set(gca,'FontSize',16);
    end
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

%% Save Results
feat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(feat_out_dir,'dir')
    mkdir(feat_out_dir);
end
stat_out_fname = [feat_out_dir SBJ_id '_' feat_id '_' proc_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','SBJs','ft_amp','pk_amp','pk_times','miss_fts','bad_pks');

end
