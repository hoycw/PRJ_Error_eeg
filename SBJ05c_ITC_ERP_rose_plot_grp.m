function SBJ05c_ITC_ERP_rose_plot_grp(SBJs,conditions,proc_id,itc_an_id,phs_id,erp_an_id,plt_id,save_fig,varargin)
%% Compute and plot ITPC matrix for group with ERP on top
% INPUTS:
%   conditions [str] - group of condition labels to segregate trials

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

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
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' erp_an_id '_vars.m'];
eval(an_vars_cmd);
erp_an = an;
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' itc_an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~an.complex; error('why run this without ITPC an_vars?'); end
if ~strcmp(an.event_type,erp_an.event_type); error('itc and erp event mismatch!'); end
if ~all(an.trial_lim_s==erp_an.trial_lim_s); error('itc and erp trial_lim_s mismatch!'); end
% phs_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' phs_id '_vars.m'];
% eval(phs_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select conditions (and trials)
[grp_lab, ~, ~, ~] = fn_group_label_styles(conditions);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(conditions);
% if ~strcmp(st.model_lab,{'DifOut','Out'}); error('not ready for surprise trials!'); end
grp_cond_lab = cell(size(grp_lab));
for grp_ix = 1:numel(grp_lab)
    [grp_cond_lab{grp_ix}, ~, ~, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
end

%% Load Peak Timing Information
if strcmp(phs_id,'FRN') && all(isfield(st,{'pk_reg_id','pk_stat_id','pk_an_id'}))
    % Load previous stats
    tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' st.pk_stat_id '_' st.pk_an_id '.mat']);
    if numel(tmp.SBJs)~=numel(SBJs) || ~all(strcmp(tmp.SBJs,SBJs))
        error(['Not same SBJs in ' stat_id ' and ' st.pk_stat_id]);
    end
    
    % Obtain peak times for target regressor
    reg_ix = find(strcmp(reg_lab,st.pk_reg_id));
    pk_ts = nan(size(tmp.time_vec));
    for t_ix = 1:numel(tmp.time_vec)
        pk_ts(t_ix) = tmp.lme{t_ix}.Coefficients.Estimate(reg_ix+1);
    end
    [~,pk_ix] = max(abs(pk_ts));
    reg_pk_time = tmp.time_vec(pk_ix);
    st.stat_lim = st.stat_lim+reg_pk_time;
else
    reg_pk_time = nan;
end

%% Load data
cfgs = [];
cfgs.avgoverfreq = 'yes';
cfgs.frequency   = phs.freq;
cfgs.latency     = phs.lim;
for s = 1:numel(SBJs)
    fprintf('-------------------- Processing %s ------------------------\n',SBJs{s});
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    load([SBJ_vars.dirs.proc SBJs{s} '_' erp_an_id '.mat'],'roi');
    load([SBJ_vars.dirs.proc SBJs{s} '_' proc_id '_' itc_an_id '.mat'],'tfr');
    if s==1
        time_vec = tfr.time;
        roi_time_vec = roi.time{1};
        fois     = tfr.freq;
        freq_idx = tfr.freq>=phs.freq(1) & tfr.freq<=phs.freq(2);
        time_idx = tfr.time>=phs.lim(1) & tfr.time<=phs.lim(2);
        ch_list  = tfr.label;
        itpc_all = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(fois) numel(time_vec)]);
        means    = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(roi_time_vec)]);
        ang      = cell([numel(cond_lab) numel(ch_list)]);
    end
    
    % Compute ITPC
    cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        % Phase angle computation
        cfgs.trials = cond_trial_ix;
        complex = ft_selectdata(cfgs,tfr);
        for ch_ix = 1:numel(ch_list)
            tmp_ang = angle(complex.fourierspctrm(:,ch_ix,:,:));
            ang{cond_ix,ch_ix} = [ang{cond_ix,ch_ix}; tmp_ang(:)];
        end
        
        % ITPC Computation
        F = tfr.fourierspctrm(cond_trial_ix,:,:,:);
        itpc = F./abs(F);                               % Normalize to unit circle
        itpc = sum(itpc,1);                             % Sum phase angles
        itpc = abs(itpc)/numel(cond_trial_ix);        % Get mean of angles for consistency
        % Add to running total
        itpc_all(cond_ix,s,:,:,:) = squeeze(itpc);
        
        % ERP Computation
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(roi_time_vec)]);
        for t_ix = 1:numel(cond_trial_ix)
            trials(:,t_ix,:) = roi.trial{cond_trial_ix(t_ix)};
        end
        means(cond_ix,s,:,:) = mean(trials,2);
    end
    clear bhv roi tfr itpc SBJ_vars cond_idx
end

% Normalize by number of SBJs
itpc_avg = nan([numel(cond_lab) numel(ch_list) numel(fois) numel(time_vec)]);
itpc_avg(:,:,:,:) = nanmean(itpc_all,2);    % squeeze will take out ch dimension

%% Get event timing for plotting
evnt_times = zeros(size(plt.evnt_lab));
if strcmp(an.event_type,'S')
    for evnt_ix = 1:numel(plt.evnt_lab)
        switch plt.evnt_lab{evnt_ix}
            case 'S'
                evnt_times(evnt_ix) = 0;
            case 'R'
                evnt_times(evnt_ix) = prdm_vars.target;
            case {'F','Fon'}
                evnt_times(evnt_ix) = prdm_vars.target+prdm_vars.fb_delay;
            case 'Foff'
                evnt_times(evnt_ix) = prdm_vars.target+prdm_vars.fb_delay+prdm_vars.fb;
            otherwise
                error(['Unknown event type in plt: ' plt.evnt_lab{evnt_ix}]);
        end
    end
elseif strcmp(an.event_type,'F')
    evnt_times(1) = 0;
else
    error('Unknown an.event_type');
end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' itc_an_id '/' conditions '/' erp_an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Create a figure for each channel
for ch_ix = 1:numel(ch_list)
    %% Compute plotting data    
    % Compute means and variance
    itpc_phase_win = nan(size(cond_lab));
    plt_means = NaN([numel(cond_lab) numel(roi_time_vec)]);
    sems      = NaN([numel(cond_lab) numel(roi_time_vec)]);
    for cond_ix = 1:numel(cond_lab)
        plt_means(cond_ix,:) = squeeze(mean(means(cond_ix,:,ch_ix,:),2));
        sems(cond_ix,:) = squeeze(std(means(cond_ix,:,ch_ix,:),[],2))./sqrt(numel(SBJs))';
        
        itpc_phase_win(cond_ix) = squeeze(mean(mean(itpc_avg(cond_ix,ch_ix,freq_idx,time_idx),4),3));
    end
    
    %% Create plot
    fig_name = ['GRP_' conditions '_' itc_an_id '_' phs_id '_' erp_an_id '_' ch_list{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Get color lims per condition
    clim = zeros([numel(cond_lab) 2]);
    for cond_ix = 1:numel(cond_lab)
        vals = itpc_avg(cond_ix,ch_ix,:,:);
        clim(cond_ix,:) = [min(vals(:)) max(vals(:))];
    end
    tick_ix = 1:3:numel(fois);
    yticklab = cell(size(tick_ix));
    for f = 1:numel(tick_ix)
        yticklab{f} = num2str(fois(tick_ix(f)),'%.1f');
    end
    
    % Condition Plots
    for cond_ix = 1:length(cond_lab)
%         subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix);
        subplot(2,numel(cond_lab),cond_ix);
        % Plot ITC Matrix
        yyaxis left
        %contourf(time_vec, fois, squeeze(itpc_avg(cond_ix,ch_ix,:,:)));
        imagesc(time_vec, fois, squeeze(itpc_avg(cond_ix,ch_ix,:,:)));% 1:numel(fois)
        set(gca,'YDir','normal');
%         set(gca,'YTick',1:3:numel(fois));
%         set(gca,'YTickLabels',yticklab);
        ylabel('Frequency (Hz)');
        caxis([min(clim(:,1)) max(clim(:,2))]);
        colorbar('northoutside');
        
        % Plot time-frequency ROI
        line([phs.lim(1) phs.lim(1)], phs.freq, 'Color','k','LineWidth',2,'LineStyle','--');
        line([phs.lim(2) phs.lim(2)], phs.freq, 'Color','k','LineWidth',2,'LineStyle','--');
        line(phs.lim, [phs.freq(1) phs.freq(1)], 'Color','k','LineWidth',2,'LineStyle','--');
        line(phs.lim, [phs.freq(2) phs.freq(2)], 'Color','k','LineWidth',2,'LineStyle','--');
        
        % Plot Events
        for evnt_ix = 1:numel(plt.evnt_lab)
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        % Plot Means (and variance)
        yyaxis right
        shadedErrorBar(roi_time_vec, plt_means(cond_ix,:), sems(cond_ix,:),...
                'lineProps',{'Color','k','LineWidth',2,...
                'LineStyle','-'},'patchSaturation',0.3);
        ylabel('Amplitude (uV)');
        
        % Axis Parameters
        title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        set(gca,'FontSize',16);
        axis tight;
        
        % Plot Polar Historgram
        subplot(2,numel(cond_lab),cond_ix+numel(cond_lab));
        polarhistogram(ang{cond_ix,ch_ix},[-pi:pi/5:pi],'Normalization','probability');
        title([cond_lab{cond_ix} ': ITPC = ' num2str(itpc_phase_win(cond_ix))]);
    end
    
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.png'];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end

end
