function SBJ06a_OB_ERP_save_mean_window(SBJ_id,proc_id,feat_id,varargin)
%% Save amplitude and latency of condition-averaged ERP peaks or mean window in oddball (OB) task
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load and average single-trial ERP within condition
%   Find feature based on peaks, take largest amplitude/latency
%       NOTE: feature measure is currently constrained to the same across features
%   Plot ERPs and detected peaks
%   Save amplitudes and latencies
% INPUTS:
%   SBJ_id [str]  - ID of subject list for group
%   proc_id [str] - ID of oddball preprocessing pipeline
%   feat_id [str] - ID of the feature extraction parameters to use
%       ft.an_id  = 'ERP_all_S2t1_dm2t0_fl05t20'
%       ft.grp_id = conditions to extract features (likely 'Tar','Odd','rare')
%       ft.measure:
%           'grpMW'- Mean window amplitude based on group peak and latency
%           'sbjMW'- Mean window amplitude based on single SBJ peak and latency
%           'sbjPk'- Single SBJ peaks and latencies
%   varargin:
%       plot_feat [0/1] - binary flag for plotting each ERP with peaks overlay
%           default: 1
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)
%   ft_amp [float] - amplitudes of [positive, negative] features
%   ft_times [float] - times (in sec) of [positive, negative] features

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

%% Load Data 
if ~contains(proc_id,'odd'); error('proc_id must be for oddball task!'); end
feat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/feat_vars/' feat_id '_vars.m'];
eval(feat_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(ft.grp_id);

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
        grp_erps = nan([numel(cond_lab) numel(ch_list) numel(time_vec)]);
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

% Compute group (grand-average) ERPs while keeping condition/channel dims
if strcmp(ft.measure,'grpMW')
    for cond_ix = 1:numel(cond_lab)
        grp_erps(cond_ix,:,:) = mean(erps(cond_ix,:,:,:),2);
    end
end

%% Compute Peak Amplitude and Latency
ft_amp   = nan([numel(SBJs) numel(ft.name)]);
ft_times = nan([numel(SBJs) numel(ft.name)]);
miss_fts = false([numel(SBJs) numel(ft.name)]);
for ft_ix = 1:numel(ft.name)
    % Select peak feature parameters
    cond_ix = find(strcmp(cond_lab,ft.cond{ft_ix}));
    ch_ix = find(strcmp(ch_list,ft.chan{ft_ix}));
    
    % Find Peak Search Range
    ft_rng = zeros([2 1]);
    for lim_ix = 1:2
        % Match to time vector
        [~, ft_rng(lim_ix)] = min(abs(time_vec-ft.lim(ft_ix,lim_ix)));
    end
    ft_time_vec = time_vec(ft_rng(1):ft_rng(2));
    if ft.pk_sign(ft_ix)~=1 && ft.pk_sign(ft_ix)~=-1
        error('ft.pk_sign not 1/-1!');
    end
    
    % Detect ERP feature and latency
    if strcmp(ft.measure,'grpMW')
        % Obtain group ERP peak time in window
        [~,pk_ix] = max(squeeze(grp_erps(cond_ix,ch_ix,ft_rng(1):ft_rng(2)))*ft.pk_sign(ft_ix));
        ft_times(:,ft_ix) = repmat(ft_time_vec(pk_ix),size(SBJs));
        
        % Compute mean window
        for s = 1:numel(SBJs)
            [~, ft_win_start] = min(abs(time_vec-(ft_times(s,ft_ix)+ft.mn_lim(ft_ix,1))));
            [~, ft_win_end]   = min(abs(time_vec-(ft_times(s,ft_ix)+ft.mn_lim(ft_ix,2))));
            ft_amp(s,ft_ix) = mean(erps(cond_ix,s,ch_ix,ft_win_start:ft_win_end));
        end
    else
        % Detect Single SBJ Peaks
        for s = 1:numel(SBJs)
            % Find all possible peak amplitudes and latencies
            [tmp_amp, tmp_lat] = findpeaks(ft.pk_sign(ft_ix)*...
                squeeze(erps(cond_ix,s,ch_ix,ft_rng(1):ft_rng(2))));
            % Convert latencies to time
            if ~isempty(tmp_lat)
                tmp_times = ft_time_vec(tmp_lat);
            end
            
            % Check for missing peaks
            if isempty(tmp_amp)
                miss_fts(s,ft_ix) = true;
                fprintf(2,'\tNo %s peak detected in %s %s!\n',...
                    ft.name{ft_ix},SBJs{s},cond_lab{cond_ix});
            else
                % Select largest peak
                [~, max_ix] = max(tmp_amp);
                ft_times(s,ft_ix) = tmp_times(max_ix);
                if strcmp(ft.measure,'sbjPk')
                    ft_amp(s,ft_ix) = ft.pk_sign(ft_ix)*tmp_amp(max_ix); % flip sign back if needed
                elseif strcmp(ft.measure,'sbjMW')
                    [~, ft_win_start] = min(abs(time_vec-(ft_times(s,ft_ix)+ft.mn_lim(ft_ix,1))));
                    [~, ft_win_end]   = min(abs(time_vec-(ft_times(s,ft_ix)+ft.mn_lim(ft_ix,2))));
                    ft_amp(s,ft_ix) = mean(erps(cond_ix,s,ch_ix,ft_win_start:ft_win_end));
                else
                    error(['Unknown ft.measure = ' ft.measure]);
                end
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
            
            if strcmp(ft.measure,'sbjPk')
                % Plot Peaks
                scatter(ft_times(s,ft_ix),ft_amp(s,ft_ix), 'o', ...
                    'MarkerEdgeColor', sbj_colors(s,:));
            else
                % Plot mean window as a line
                line(ft.mn_lim(ft_ix,:)+ft_times(s,ft_ix),repmat(ft_amp(s,ft_ix),[1 2]),...
                    'Color',sbj_colors(s,:),'LineWidth',2);
            end
        end
        
        % Plot parameters
        title([ft.name{ft_ix} ': ' ft.cond{ft_ix} ', ' ft.chan{ft_ix}]);
        xlabel('Time (s)'); ylabel('Amplitude (uV)');
        xlim([ft.plot_lim(1) ft.plot_lim(2)]);
        set(gca,'FontSize',16);
    end
    
    % Save figure
    if save_fig
        fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' ft.an_id '/OB_feat_ts/'];
        if ~exist(fig_dir,'dir') && save_fig
            mkdir(fig_dir);
        end
        
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
save(stat_out_fname,'-v7.3','SBJs','ft_amp','ft_times','miss_fts');

end
