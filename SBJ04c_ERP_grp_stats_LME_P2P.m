function SBJ04c_ERP_grp_stats_LME_P2P(SBJ_id,proc_id,an_id,stat_id,varargin)
%% Compute linear mixed-effects model on peak-to-peak FRN for condition-averaged ERPs
%   Only runs for one channel
%   Early versions used jackknife to deal with waveshape variability, but
%       abandoned because DF correction for F stat was unclear, and missing
%       FRN needed to be zero instead of NaN (couldn't correct with missing data)
% COMPUTATIONS:
%   Select trials for conditions of interest
%   Load and average single-trial ERP and design matrix (model regressors, SBJ factor) within condition
%       Optional: z-score model regressors within SBJ
%   Compute and plot group concatenated model (design matrix) and correlations
%   Find positive and negative peaks:
%       Must find positive and negative peaks in their respective windows
%       Positive peak must precede negative peak
%   Compute amplitude differences:
%       For multiple peak pairs, selects largest amplitude difference
%       Rejects pairs with amplitude difference in wrong direction (negative > positive)
%   Plot peak detection errors
%       Optional: Plot ERPs and detected peaks
%   Run linear mixed effects model per time point or per electrode
%   Correct for multiple comparisons (FDR for regressors)
% INPUTS:
%   SBJ_id [str] - ID of subject list for group
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
%   varargin:
%       plot_erps [0/1] - binary flag for plotting each ERP with peaks overlay
%           default: 0
%       plot_peaks [0/1] - binary flag for peak detection errors and summary plot of both peaks
%           default: 1
%       fig_vis [str] - {'on','off'} to visualize figure on desktop
%           default: 'on'
%       fig_ftype [str] - file extension for saving fig
%           default: 'png'
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class, one cell per channel
%   qvals [float array] - [n_regressors, n_chan/n_time] p values adjusted for multiple comparisons 
%   SBJs [cell array] - list of SBJs used in this analysis (for double checks)
%   time_vec [float array] - time points for each model run (same length as lme)
%   ch_list [cell array] - list of channels in analysis (for double checks)
%   pk_amp [float] - amplitudes of [positive, negative] peaks
%   pk_times [float] - times (in sec) of [positive, negative] peaks
%   data [floats] - peak-to-peak amplitude differences per [condition, subject, channel]
%       contains NaNs for any SBJ condition if any peak-to-peak criteria not met

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
        elseif strcmp(varargin{v},'plot_erps')
            plot_erps = varargin{v+1};
        elseif strcmp(varargin{v},'plot_peaks')
            plot_peaks = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');    fig_vis = 'on'; end
if ~exist('fig_ftype','var');  fig_ftype = 'png'; end
if ~exist('save_fig','var');   save_fig = 1; end
if ~exist('plot_erps','var');  plot_erps = 0; end
if ~exist('plot_peaks','var'); plot_peaks = 1; end
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' stat_id '/GLM_inputs/'];
if ~exist(fig_dir,'dir') && save_fig
    mkdir(fig_dir);
end

%% Load Data 
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.measure,'p2p') || ~strcmp(st.an_style,'lme'); error('This script should only be run on LME+p2p!'); end
if plot_peaks && ~strcmp(st.measure,'p2p'); error('why plot peaks if not using p2p?'); end
if strcmp(st.grp_method,'jackknife'); error('jackknife not figured out for GLM/LME!'); end
if ~strcmp(st.grp_method,'sbj'); error(['Unknown st.grp_method: ' st.grp_method]); end

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Get model and condition parameters
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, cond_colors, ~, ~] = fn_condition_label_styles(st.stat_cond);

%% Load Behavior
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

%% Load Data and Build Model
cfgs  = []; cfgs.latency = st.stat_lim;
model = zeros([numel(cond_lab)*numel(SBJs) numel(reg_lab)]);
sbj_factor = [];
full_model_ix = 0;
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',an_id,'.mat'],'roi');
    if numel(roi.label)>1; error('only ready for one channel in P2P GLM!'); end
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_roi = ft_selectdata(cfgs, roi);
    
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = st_roi.time{1};
        ch_list  = st_roi.label;
        erps  = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
    end
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' st.model_id '.mat']);
    
    % Z-score SBJ model regressors
    sbj_model = NaN([sum(full_cond_idx{s}~=0) size(tmp.model,2)]);
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            sbj_model(:,reg_ix) = ...
                (tmp.model(full_cond_idx{s}~=0,reg_ix)-nanmean(tmp.model(full_cond_idx{s}~=0,reg_ix)))./...
                nanstd(tmp.model(full_cond_idx{s}~=0,reg_ix));
        end
    else
        sbj_model = tmp.model(full_cond_idx{s}~=0,:);
    end
    
    % Compute ERPs and average model within condition
    cond_idx = fn_condition_index(cond_lab, bhvs{s});
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        % Compute ERP within condition
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for trl_ix = 1:numel(cond_trial_ix)
            trials(:,trl_ix,:) = st_roi.trial{cond_trial_ix(trl_ix)};
        end
        erps(cond_ix,s,:,:) = mean(trials,2);
        
        % Average model within condition
        full_model_ix = full_model_ix + 1;
        model(full_model_ix, :) = mean(sbj_model(cond_trial_ix,:),1);
    end
    
    % Track SBJ in design matrix
    sbj_factor = [sbj_factor; s*ones([numel(cond_lab) 1])];
    
    clear tmp roi st_roi sbj_model trials cond_idx cond_trial_ix
end

%% Compute and plot correlations between regressors
reg_corr = corr(model,'rows','complete');
rvals = reshape(triu(reg_corr,1),[numel(reg_corr) 1]);
[~,max_r_ix] = max(abs(rvals));
if st.z_reg
    reg_vifs = fn_variance_inflation_factor(model);
else
    reg_vifs = fn_variance_inflation_factor(zscore(model));
end
[max_vif,max_vif_ix] = max(reg_vifs);

% Create figure directory
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
fig_dir = [stat_out_dir st.model_id '_' st.stat_cond '_erpP2P_plots/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Plot design matrix
fig_name = [SBJ_id '_' st.model_id '_' st.stat_cond '_design'];
figure('Name',fig_name);
imagesc(model);
xticklabels(reg_lab);
colorbar;
saveas(gcf,[fig_dir fig_name '.png']);

% Plot regressor correlation matrix
fig_name = [SBJ_id '_' st.model_id '_' st.stat_cond '_design_corr_VIFs'];
figure('Name',fig_name);
imagesc(reg_corr);
set(gca,'XLim',[0.5 numel(reg_lab)+0.5]);
set(gca,'XTick',1:numel(reg_lab));
set(gca,'XTickLabels',reg_lab);
set(gca,'YLim',[0.5 numel(reg_lab)+0.5]);
set(gca,'YTick',1:numel(reg_lab));
set(gca,'YTickLabels',reg_lab);
colorbar;
title(['Max r = ' num2str(rvals(max_r_ix),'%.3f')]);
set(gca,'FontSize',16);

% Plot VIFs
subplot(1,2,2);
bars = bar(1:numel(reg_lab),reg_vifs);
bars.FaceColor = 'flat';
bars.CData(2,:) = [.5 0 .5];
for reg_ix = 1:numel(reg_lab)
    if reg_vifs(reg_ix) >= 5 && reg_vifs(reg_ix) < 10
        bars.CData(reg_ix,:) = [1 .5 0];
    elseif reg_vifs(reg_ix) >= 10
        bars.CData(reg_ix,:) = [1 0 0];
    else
        bars.CData(reg_ix,:) = [0 0 0];
    end
end
set(gca,'XLim',[0.5 numel(reg_lab)+0.5]);
set(gca,'XTick',1:numel(reg_lab));
set(gca,'XTickLabels',reg_lab);
set(gca,'FontSize',16);
title(['Max VIF: ' reg_lab{max_vif_ix} '=' num2str(max_vif,'%.2f')]);
saveas(gcf,[fig_dir fig_name '.png']);

%% Compute Peak-to-Peak Data
% Find Peak Search Range
pk_rng = zeros([2 2]);
pk_sign_str = cell([2 1]);
% If any stat_id changes second peak from negative, check code before running!
if st.pk_sign(2)~=-1; error('second peak is not negative!'); end
for pk_ix = 1:2
    for lim_ix = 1:2
        % Match to time vector
        [~, pk_rng(pk_ix,lim_ix)] = min(abs(time_vec-st.pk_lim(pk_ix,lim_ix)));
    end
    if st.pk_sign(pk_ix)==1
        pk_sign_str{pk_ix} = 'positive';
    elseif st.pk_sign(pk_ix)==-1
        pk_sign_str{pk_ix} = 'negative';
    else
        error('st.pk_sign not 1/-1!');
    end
end

% Detect peaks and compute amplitude difference
data     = nan([numel(cond_lab) numel(SBJs) numel(ch_list)]);
pk_amp   = nan([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
pk_times = nan([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
miss_pks = false([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
bad_pks  = false([numel(cond_lab) numel(SBJs) numel(ch_list)]);
for s = 1:numel(SBJs)
    for cond_ix = 1:numel(cond_lab)
        for ch_ix = 1:numel(ch_list)
            % Find all possible positive and negative peak amplitudes and latencies
            tmp_amp = cell([2 1]); tmp_times = cell([2 1]);
            for pk_ix = 1:2
                [tmp_amp{pk_ix}, tmp_lat] = findpeaks(st.pk_sign(pk_ix)*...
                    squeeze(erps(cond_ix,s,ch_ix,pk_rng(pk_ix,1):pk_rng(pk_ix,2))));
                % Convert latencies to time
                if ~isempty(tmp_lat)
                    tmp_times{pk_ix} = time_vec(tmp_lat+pk_rng(pk_ix,1)-1);
                end
            end
            
            % Check for missing peaks
            if isempty(tmp_amp{1}) || isempty(tmp_amp{2})
                miss_pks(cond_ix,s,ch_ix) = true;
                if isempty(tmp_amp{1})
                    fprintf(2,'\tNo pos peak detected for %s %s!\n',...
                        SBJs{s},cond_lab{cond_ix});
                end
                if isempty(tmp_amp{2})
                    fprintf(2,'\tNo neg peak detected for %s %s!\n',...
                        SBJs{s},cond_lab{cond_ix});
                end
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
                        % !!! should pk_ix here be hard coded as 2???
                        if tmp_amp{1}(pre_neg_pos_ix(i))+tmp_amp{2}(neg_pk_ix) > tmp_p2p(neg_pk_ix)
                            tmp_p2p(neg_pk_ix) = tmp_amp{1}(pre_neg_pos_ix(i))+tmp_amp{2}(neg_pk_ix);
                            p2p_pos_ix(neg_pk_ix) = pre_neg_pos_ix(i);
                        end
                    end
                end
                
                % Select pair of peaks with max peak-to-peak difference
                [~, max_p2p_neg_ix] = max(tmp_p2p);
                pk_amp(cond_ix,s,ch_ix,1)   = tmp_amp{1}(p2p_pos_ix(max_p2p_neg_ix));
                pk_times(cond_ix,s,ch_ix,1) = tmp_times{1}(p2p_pos_ix(max_p2p_neg_ix));
                pk_amp(cond_ix,s,ch_ix,2)   = -tmp_amp{2}(max_p2p_neg_ix);  % Flip sign back after findpeaks
                pk_times(cond_ix,s,ch_ix,2) = tmp_times{2}(max_p2p_neg_ix);
            end
            
            % Compute peak-to-peak amplitude difference
            if ~miss_pks(cond_ix,s,ch_ix)
                amp_diff = pk_amp(cond_ix,s,ch_ix,1)-pk_amp(cond_ix,s,ch_ix,2);
                % Toss differences that go in the wrong direction or order
                %   This indicates a complex, multi-peak waveform...
                if amp_diff > 0
                    data(cond_ix,s,ch_ix) = amp_diff;
                else
                    bad_pks(cond_ix,s,ch_ix) = true;
                    fprintf(2,'\tBad peaks detected for %s %s!\n',...
                        SBJs{s},cond_lab{cond_ix});
                end
            end
            
            % Plot Peak Detection Errors (missing or bad peaks)
            if plot_peaks && (any(miss_pks(cond_ix,s,ch_ix,:)) || bad_pks(cond_ix,s,ch_ix))
                fig_name = ['bad_peaks_' st.grp_method '_' SBJs{s} ...
                    '_' cond_lab{cond_ix} '_' ch_list{ch_ix}];
                figure('Name',fig_name,'Visible',fig_vis);
                
                % Plot ERP
                plot(time_vec(pk_rng(1,1):pk_rng(2,2)),...
                    squeeze(erps(cond_ix,s,ch_ix,pk_rng(1,1):pk_rng(2,2))),...
                    'Color',cond_colors{cond_ix});
                hold on;
                
                % Plot Peaks on top
                for pk_ix = 1:2
                    [tmp_amp, tmp_lat] = ...
                        findpeaks(st.pk_sign(pk_ix)*...
                        squeeze(erps(cond_ix,s,ch_ix,pk_rng(pk_ix,1):pk_rng(pk_ix,2))));
                    if st.pk_sign(pk_ix)==1; mrkr_color = 'b'; else mrkr_color = 'r'; end
                    if ~miss_pks(cond_ix,s,ch_ix,pk_ix)
                        scatter(time_vec(tmp_lat+pk_rng(pk_ix,1)-1),...
                            st.pk_sign(pk_ix)*tmp_amp, 'o', mrkr_color);
                    end
                end
                
                % Create error message for plot title
                bad_str = [];
                if bad_pks(cond_ix,s,ch_ix)
                    bad_str = [bad_str ' bad;'];
                end
                if any(miss_pks(cond_ix,s,ch_ix,:))
                    bad_str = [bad_str ' miss' num2str(sum(miss_pks(cond_ix,s,ch_ix,:))) ';'];
                else
                    % Plot line connecting peaks
                    line([pk_times(cond_ix,s,ch_ix,1) pk_times(cond_ix,s,ch_ix,2)],...
                        [pk_amp(cond_ix,s,ch_ix,1) pk_amp(cond_ix,s,ch_ix,2)],...
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
end

%% Plot Data
% Plot ERPs (all overlapping, colored by condition)
if plot_erps
    for ch_ix = 1:numel(ch_list)
        fig_name = ['ERPs_' st.grp_method '_' ch_list{ch_ix} '_' SBJ_id];
        figure('Name',fig_name,'Visible',fig_vis);
        
        % Get axis limits
        tmp = erps(:,:,ch_ix,:);
        ylims = [min(tmp(:)) max(tmp(:))];
        
        for cond_ix = 1:numel(cond_lab)
%             subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix); hold on;
            subplot(2,3,cond_ix); hold on;
            
            % Plot ERPs
            plot(time_vec,squeeze(erps(cond_ix,:,ch_ix,:)),...
                'Color',cond_colors{cond_ix});%,'LineStyle',cond_styles{cond_ix});
            
            % Plot parameters
            title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
            ylim(ylims);
            fn_min_white_space(gca);
        end
        
        % Save figure
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
    end
end

% Plot peaks in peak range
if plot_peaks
    for ch_ix = 1:numel(ch_list)
        fig_name = ['ERPs_peaks_' st.grp_method '_'  ch_list{ch_ix} '_' SBJ_id];
        figure('Name',fig_name,'Visible',fig_vis);
        
        % Get axis limits
        tmp = pk_amp(:,:,ch_ix,:);
        ylims = [min(tmp(:)) max(tmp(:))];
        tmp = pk_times(:,:,ch_ix,:);
        xlims = [min(tmp(:)) max(tmp(:))];
        
        for cond_ix = 1:numel(cond_lab)
%             subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix); hold on;
            subplot(2,3,cond_ix); hold on;
            
            % Plot ERPs
            plot(time_vec(pk_rng(1,1):pk_rng(2,2)),...
                squeeze(erps(cond_ix,:,ch_ix,pk_rng(1,1):pk_rng(2,2))),...
                'Color',cond_colors{cond_ix});%,'LineStyle',cond_styles{cond_ix});
            
            % Plot Peaks
            scatter(pk_times(cond_ix,:,ch_ix,1), pk_amp(cond_ix,:,ch_ix,1), 'o', 'b');
            scatter(pk_times(cond_ix,:,ch_ix,2), pk_amp(cond_ix,:,ch_ix,2), 'o', 'r');
            
            % Plot line connecting peaks
            for s = 1:numel(SBJs)
                line([pk_times(cond_ix,s,ch_ix,1) pk_times(cond_ix,s,ch_ix,2)],...
                    [pk_amp(cond_ix,s,ch_ix,1) pk_amp(cond_ix,s,ch_ix,2)],...
                    'Color','k');
            end
            
            % Plot parameters
            title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
            ylim(ylims);
            xlim(xlims);
            fn_min_white_space(gca);
        end
        
        % Save figure
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
    end
end

%% Compute Statistics
fprintf('========================== Running Stats ==========================\n');
tic

% Build Model Table
tbl = table;
for reg_ix = 1:numel(reg_lab)
    tbl.(reg_lab{reg_ix}) = model(:,reg_ix);
end
tbl.SBJ = categorical(sbj_factor);

% Create Model Formula
reg_formula = strjoin(reg_lab,' + ');
formula = ['ERP ~ ' reg_formula ' + (1|SBJ)'];  % random intercepts for SBJ

% Reshape means for LME input
st_data = reshape(data,[numel(cond_lab)*numel(SBJs), numel(ch_list)]);
if any(isnan(st_data(:))); warning([num2str(sum(isnan(st_data(:)))) ' NaNs in LME data!']);end

% Run Model
lme   = cell(size(ch_list));
pvals = nan([numel(reg_lab) numel(ch_list)]);
for ch_ix = 1:numel(ch_list)
    tbl.ERP    = st_data(:,ch_ix);
    lme{ch_ix} = fitlme(tbl,formula);
    pvals(:,ch_ix) = lme{ch_ix}.Coefficients.pValue(2:end); % Skip intercept
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

fprintf('\t\t Stats Complete:');
toc

%% Save Results
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(stat_out_dir,'dir')
    mkdir(stat_out_dir);
end
stat_out_fname = [stat_out_dir SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals','SBJs','time_vec','ch_list','pk_amp','pk_times','data');

end
