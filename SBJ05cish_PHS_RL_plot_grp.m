function SBJ05cish_PHS_RL_plot_grp(SBJs,proc_id,itc_an_id,phs_id,stat_id,plt_id,save_fig,varargin)
%% Compute and plot ITPC matrix for group with ERP on top

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
addpath([app_dir 'CircStat/']);
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

%% Load Parameters
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' erp_an_id '_vars.m'];
eval(an_vars_cmd);
erp_an = an;
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' itc_an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~an.complex; error('why run this without ITPC an_vars?'); end
if ~strcmp(an.event_type,erp_an.event_type); error('itc and erp event mismatch!'); end
if ~all(an.trial_lim_s==erp_an.trial_lim_s); error('itc and erp trial_lim_s mismatch!'); end
phs_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' phs_id '_vars.m'];
eval(phs_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select conditions (and trials)
model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[grp_lab, ~, ~, ~] = fn_group_label_styles(conditions);
[cond_lab, ~, ~, ~, ~] = fn_condition_label_styles(conditions);
% if ~strcmp(st.model_lab,{'DifOut','Out'}); error('not ready for surprise trials!'); end
grp_cond_lab = cell(size(grp_lab));
for grp_ix = 1:numel(grp_lab)
    [grp_cond_lab{grp_ix}, ~, ~, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
end

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

%% Load Peak Timing Information
% if strcmp(st.measure,'mean') && all(isfield(st,{'pk_reg_id','pk_stat_id','pk_an_id'}))
%     % Load previous stats
%     tmp = load([root_dir 'PRJ_Error_eeg/data/GRP/GRP_' st.pk_stat_id '_' st.pk_an_id '.mat']);
%     if numel(tmp.SBJs)~=numel(SBJs) || ~all(strcmp(tmp.SBJs,SBJs))
%         error(['Not same SBJs in ' stat_id ' and ' st.pk_stat_id]);
%     end
%     
%     % Obtain peak times for target regressor
%     reg_ix = find(strcmp(reg_lab,st.pk_reg_id));
%     pk_ts = nan(size(tmp.time_vec));
%     for t_ix = 1:numel(tmp.time_vec)
%         pk_ts(t_ix) = tmp.lme{t_ix}.Coefficients.Estimate(reg_ix+1);
%     end
%     [~,pk_ix] = max(abs(pk_ts));
%     reg_pk_time = tmp.time_vec(pk_ix);
%     st.stat_lim = st.stat_lim+reg_pk_time;
% end

%% Load data
cfgs = [];
cfgs.avgoverfreq = 'yes';
cfgs.frequency   = phs.freq;
cfgs.latency     = phs.lim;

model = zeros([sum(n_trials) numel(reg_lab)]);
sbj_factor  = zeros([sum(n_trials) 1]);
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
        ang      = cell(size(ch_list));
    
        % Load RL Model
        tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
        % model(1:n_trials(s),:) = tmp.model;
        % Z-score SBJ model regressors
        sbj_model = NaN(size(tmp.model));
        if st.z_reg
            for reg_ix = 1:numel(reg_lab)
                sbj_model(:,reg_ix) = ...
                    (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
            end
        else
            sbj_model = model;
        end
        model(1:n_trials(s),:) = sbj_model;
        
        % Track SBJ
        sbj_factor(1:n_trials(s),end) = s*ones([n_trials(s) 1]);
    else
        % Load RL Model
        tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
        % model(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),:) = tmp.model;
        % Z-score SBJ model regressors
        sbj_model = NaN(size(tmp.model));
        if st.z_reg
            for reg_ix = 1:numel(reg_lab)
                sbj_model(:,reg_ix) = ...
                    (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
            end
        else
            sbj_model = model;
        end
        model(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),:) = sbj_model;
        
        % Track SBJ
        sbj_factor(sum(n_trials(1:s-1))+1:sum(n_trials(1:s))) = s*ones([n_trials(s) 1]);
    end
    
    % Compute Mean Phase Angle per trial
    complex = ft_selectdata(cfgs,tfr);
    for ch_ix = 1:numel(ch_list)
        ang{ch_ix} = [ang{ch_ix}; squeeze(angle(nanmean(complex.fourierspctrm(:,ch_ix,:,:),4)))];
    end
        
    clear bhv roi tfr itpc SBJ_vars cond_idx
end

% Normalize by number of SBJs
% itpc_avg = nan([numel(cond_lab) numel(ch_list) numel(fois) numel(time_vec)]);
% itpc_avg(:,:,:,:) = nanmean(itpc_all,2);    % squeeze will take out ch dimension

% % Concat SBJs
% for ch_ix = 1:numel(ch_list)
%     phs = vertcat(ang{:,ch_ix});
% end

%% Phase-Model Correlations
phs_corr = nan(size(reg_lab));
phs_pval = nan(size(reg_lab));
for reg_ix = 1:numel(reg_lab)
    for ch_ix = 1:numel(ch_list)
        [phs_corr(reg_ix), phs_pval(reg_ix)] = circ_corrcl(ang{ch_ix}, model(:,reg_ix));
    end
end

%% Bin regressor by phase
bins = [-pi:pi/4:pi];
reg_bar = nan([numel(reg_lab) numel(bins)]);
for reg_ix = 1:numel(reg_lab)
    for b_ix = 2:numel(bins)
        phs_idx = ang{ch_ix}>=bins(b_ix-1) & ang{ch_ix}<bins(b_ix);
        reg_bar(reg_ix,b_ix) = mean(model(phs_idx,reg_ix));
    end
end

% plot
figure;
for r = 1:7
    subplot(4,2,r);
    bar(bins,reg_bar(r,:));
    title([reg_lab{r} ': r = ' num2str(phs_corr(r)) ' (p = ' num2str(phs_pval(r)) ')']);
end

%% Or rose plot by binned regressor?
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
