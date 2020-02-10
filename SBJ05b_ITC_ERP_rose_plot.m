function SBJ05b_ITC_ERP_rose_plot(SBJ,conditions,proc_id,itc_an_id,erp_an_id,plt_id,save_fig,varargin)
%% Compute and plot ITPC matrix for single SBJ; plot ERP on top
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
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' erp_an_id '_vars.m'];
eval(an_vars_cmd);
erp_an = an;
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' itc_an_id '_vars.m'];
eval(an_vars_cmd);
if an.avgoverfreq; error('why run this with only 1 freq in an_vars?'); end
if ~an.itpc; error('why run this without ITPC an_vars?'); end
if ~strcmp(an.event_type,erp_an.event_type); error('itc and erp event mismatch!'); end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Load data
load([SBJ_vars.dirs.proc SBJ '_' erp_an_id '.mat']);
load([SBJ_vars.dirs.proc SBJ '_' proc_id '_' itc_an_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

% Select conditions (and trials)
[grp_lab, ~, ~] = fn_group_label_styles(conditions);
[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);
% if ~strcmp(st.model_lab,{'DifOut','Out'}); error('not ready for surprise trials!'); end
grp_cond_lab = cell(size(grp_lab));
for grp_ix = 1:numel(grp_lab)
    [grp_cond_lab{grp_ix}, ~, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
end
cond_idx = fn_condition_index(cond_lab, bhv);

%% Angle Extraction
an.lim = [0.18 0.3];
an.freq = [4 8];
figure;
ang = cell(size(cond_lab));
erps = cell(size(cond_lab));
trials = cell(size(cond_lab));
cfgerp = []; cfgerp.latency = an.lim;
cfgs = [];
cfgs.avgoverfreq = 'yes';
cfgs.frequency = an.freq;
cfgs.latency = an.lim;
for cond_ix = 1:numel(cond_lab)
    cfgs.trials = find(cond_idx==cond_ix);
    cfgerp.trials = find(cond_idx==cond_ix);
    complex = ft_selectdata(cfgs,tfr);
    ang{cond_ix} = squeeze(angle(complex.fourierspctrm));
    
    
    erps{cond_ix} = ft_selectdata(cfgerp,roi);
    trials{cond_ix} = nan([numel(roi.label) numel(erps{cond_ix}.trial) numel(erps{cond_ix}.time{1})]);
    for t_ix = 1:numel(erps{cond_ix}.trial)
        trials{cond_ix}(:,t_ix,:) = erps{cond_ix}.trial{t_ix};
    end
%     F = tfr.fourierspctrm(cond_idx==cond_ix,:,:,:);
%     itpc{cond_ix} = F./abs(F);       % Normalize to unit circle
%     itpc{cond_ix} = sum(itpc{cond_ix},1);     % Sum phase angles
%     itpc{cond_ix} = abs(itpc{cond_ix})/sum(cond_idx==cond_ix);     % Get mean of angles for consistency

    subplot(2,6,cond_ix);
    shadedErrorBar(erps{cond_ix}.time{1},squeeze(mean(trials{cond_ix}(1,:,:),2)),...
        squeeze(std(trials{cond_ix}(1,:,:),[],2))./sqrt(size(trials{cond_ix},2)));
    xlim(an.lim);
    title([cond_lab{cond_ix} ': ' num2str(size(trials{cond_ix},2))]);
    
    subplot(2,6,cond_ix+6);
    polarhistogram(ang{cond_ix},[-pi:pi/5:pi],'Normalization','probability');
    title(cond_lab{cond_ix});
end


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
for ch_ix = 1:numel(tfr.label)
    %% Compute plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(roi.time{1})]);
    sems  = NaN([numel(cond_lab) numel(roi.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(ch_ix,:,:),2));
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(ch_ix,:,:),[],2))./sqrt(size(trials{cond_ix},2))';
    end
    
    %% Create plot
    fig_name = [SBJ '_' conditions '_' itc_an_id '_' erp_an_id '_' tfr.label{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    
    % Get color lims per condition
    clim = zeros([numel(cond_lab) 2]);
    for cond_ix = 1:numel(cond_lab)
        vals = itpc{cond_ix}(1,ch_ix,:,:);
        clim(cond_ix,:) = [min(vals(:)) max(vals(:))];
    end
    tick_ix = 1:3:numel(tfr.freq);
    yticklab = cell(size(tick_ix));
    for f = 1:numel(tick_ix)
        yticklab{f} = num2str(tfr.freq(tick_ix(f)),'%.1f');
    end
    
    % Condition Plots
    for cond_ix = 1:length(cond_lab)
        subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix);
        % Plot ITC Matrix
        yyaxis left
        %contourf(tfr.time, tfr.freq, squeeze(itpc{cond_ix}(1,ch_ix,:,:)));
        imagesc(tfr.time, tfr.freq, squeeze(itpc{cond_ix}(1,ch_ix,:,:)));% 1:numel(tfr.freq)
        set(gca,'YDir','normal');
%         set(gca,'YTick',1:3:numel(tfr.freq));
%         set(gca,'YTickLabels',yticklab);
        ylabel('Frequency (Hz)');
        caxis([min(clim(:,1)) max(clim(:,2))]);
        colorbar('northoutside');
        
        % Plot Events
        for evnt_ix = 1:numel(plt.evnt_lab)
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_styles{evnt_ix});
        end
        
        % Plot Means (and variance)
        yyaxis right
        shadedErrorBar(roi.time{1}, means(cond_ix,:), sems(cond_ix,:),...
                'lineProps',{'Color','k','LineWidth',2,...
                'LineStyle','-'},'patchSaturation',0.3);
        ylabel('Amplitude (uV)');
        
        title([tfr.label{ch_ix} ': ' cond_lab{cond_ix}]);
        set(gca,'XLim', [plt.plt_lim(1) plt.plt_lim(2)]);
        set(gca,'XTick', plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2));
        xlabel('Time (s)');
        set(gca,'FontSize',16);
    end
    
    % Save Figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
