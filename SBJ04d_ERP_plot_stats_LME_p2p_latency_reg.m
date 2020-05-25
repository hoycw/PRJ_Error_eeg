function SBJ04d_ERP_plot_stats_LME_p2p_latency_reg(SBJ_id,an_id,stat_id,plt_id,save_fig,varargin)
% Plots beta weights per regressor for mean window LME analyses
%   Only for single channel right now...
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
        elseif strcmp(varargin{v},'plot_data')
            plot_data = varargin{v+1};
        elseif strcmp(varargin{v},'plot_violins')
            plot_violins = varargin{v+1};
        elseif strcmp(varargin{v},'plot_latencies')
            plot_latencies = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');        fig_vis = 'on'; end
if ~exist('fig_ftype','var');      fig_ftype = 'png'; end
if ~exist('plot_data','var');      plot_data = 0; end
if ~exist('plot_violins','var');   plot_violins = 0; end
if ~exist('plot_latencies','var'); plot_latencies = 0; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Analysis and Plotting Parameters
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
if ~strcmp(st.measure,'p2p') || ~strcmp(st.an_style,'lme'); error('run only for p2p LME analyses!'); end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

% Select SBJs
SBJs = fn_load_SBJ_list(SBJ_id);

% Select Conditions of Interest
[reg_lab, reg_names, reg_colors, ~]  = fn_regressor_label_styles(st.model_lab);
[cond_lab, cond_names, cond_colors, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Stats
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat']);
if numel(ch_list)>1; error('only plotting for 1 channel in this script!'); end

%% Plot Results
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' stat_id '/' plt_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Compute plotting data
% Obtain Model Parameters
plot_betas = lme{1}.Coefficients.Estimate(2:end);
r2 = lme{1}.Rsquared.Adjusted;

% Find significant regressors
sig_reg    = false(size(reg_lab));
for reg_ix = 1:numel(reg_lab)
    if any(qvals(reg_ix,:) <= st.alpha)
        sig_reg(reg_ix) = true;
    end
end

if plot_violins
    plot_data = data;
else
    plot_means = nan(size(cond_lab));
    plot_sems  = nan(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        plot_means(cond_ix) = nanmean(data(cond_ix,:),2);
        plot_sems(cond_ix)  = nanstd(data(cond_ix,:),[],2)./sqrt(size(data,2))';
    end
end

%% Plot Peak-to-Peak Data
fig_name = [SBJ_id '_' stat_id '_' an_id '_' ch_list{1}];
if plot_data
    fig_name = [fig_name '_data'];
end
if plot_violins
    fig_name = [fig_name '_violins'];
end
figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
if plot_data
    subplot(2,1,1); hold on;
    
    if plot_violins
        violins = violinplot(squeeze(plot_data(:,:))', cond_lab, ...
            'ShowMean', true, 'ViolinAlpha', 0.3);
        
        for cond_ix = 1:numel(cond_lab)
            % Fix Mean from population estimate to sample mean
            violins(cond_ix).MeanPlot.YData = repmat(nanmean(plot_data(cond_ix,:),2), [1 2]);
            [~,mean_ix] = min(abs(violins(cond_ix).ViolinPlot.YData - nanmean(plot_data(cond_ix,:),2)));
            mean_width  = violins(cond_ix).ViolinPlot.XData(mean_ix)-cond_ix;
            violins(cond_ix).MeanPlot.XData = [cond_ix-mean_width cond_ix+mean_width];
            violins(cond_ix).MeanPlot.LineWidth = 3;
            
            % Change the colors to match condition
            %   Violin
            violins(cond_ix).ViolinColor = cond_colors{cond_ix};
            %   Box plot
            violins(cond_ix).BoxPlot.FaceColor = cond_colors{cond_ix};
            violins(cond_ix).EdgeColor = cond_colors{cond_ix};
            % 	Scatter
            %         scat_colors = zeros([numel(violins(cond_ix).ScatterPlot.XData) 3]);
            %         for s = 1:numel(SBJs)
            %             scat_colors(s,:) = cond_colors{cond_ix};
            %         end
            %         violins(cond_ix).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
            %         violins(cond_ix).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
            violins(cond_ix).ScatterPlot.CData = cond_colors{cond_ix};%scat_colors;
        end
    else
        % Plot Means as bar
        b = bar(1:numel(cond_lab),plot_means,'BarWidth',plt.bar_width,'FaceColor','flat');
        for cond_ix = 1:numel(cond_lab)
            b.CData(cond_ix,:) = cond_colors{cond_ix};
            
            % Plot SEM as line
            line([cond_ix cond_ix],[plot_means(cond_ix)-plot_sems(cond_ix) plot_means(cond_ix)+plot_sems(cond_ix)],...
                'Color','k','LineWidth',2);
        end
    end
    
    % Add label and min RT for perspective
    ax = gca;
    ax.YLabel.String = 'FRN Amplitude (uV)';
    ax.XLim          = [0 numel(cond_lab)+1];
    ax.XTick         = 1:numel(cond_lab);
    ax.XTickLabel    = cond_names;
    ax.XLabel.String = 'Conditions';
    ax.Title.String  = [ch_list{1} ' Peak-to-Peak: R2 = ' num2str(r2,'%.3f') ' (n = ' num2str(numel(SBJs)) ')'];
    set(ax,'FontSize',16');
end

%% Plot Betas
if plot_data; subplot(2,1,2); end
ax = gca; hold on;

% Find significance marker distance
sig_y = max(abs(plot_betas))*plt.sig_yfudge;

% Plot Betas
b = bar(1:numel(reg_lab),plot_betas,'BarWidth',plt.bar_width,'FaceColor','flat');
for reg_ix = 1:numel(reg_lab)
    b.CData(reg_ix,:) = reg_colors{reg_ix};
    
    % Plot Significance
    if sig_reg(reg_ix)
        if qvals(reg_ix)<=plt.sig_cut3
            x_pos = [-plt.sig_xfudge 0 plt.sig_xfudge];
        elseif qvals(reg_ix)<=plt.sig_cut2
            x_pos = [-plt.sig_xfudge plt.sig_xfudge]/2;
        else
            x_pos = 0;
        end
        y_pos = plot_betas(reg_ix)+sig_y*sign(plot_betas(reg_ix));
        scatter(x_pos+reg_ix,repmat(y_pos,size(x_pos)),...
            plt.sig_sz,plt.sig_color,plt.sig_mrkr);
    end
end

% Axes and Labels
ax.YLabel.String = 'Beta Weight';
ax.XLim          = [0 numel(reg_lab)+1];
ax.XTick         = 1:numel(reg_lab);
ax.XTickLabel    = reg_names;
ax.XLabel.String = 'Regressors';
reg_str = cell(size(reg_lab));
for reg_ix = 1:numel(reg_lab)
    reg_str{reg_ix} = [reg_lab{reg_ix} ' (q=' num2str(qvals(reg_ix),'%.3f') ')'];
end
if plot_data
    ax.Title.String = strjoin(reg_str,', ');
else
    ax.Title.String = [ch_list{1} ' Peak-to-Peak: R2 = ' num2str(r2,'%.3f') ', ' strjoin(reg_str,', ')];
end
if plt.legend
    legend(b,reg_lab,'Location',plt.legend_loc);
end
ylims = ylim;
set(gca,'FontSize',16);
ax.YLim = ylims;

%% Save figure
if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    % Ensure vector graphics if saving
    if any(strcmp(fig_ftype,{'svg','eps'}))
        set(gcf, 'Renderer', 'painters');
    end
    saveas(gcf,fig_fname);
end

%% Plot Peak Latencies
if plot_latencies
    fig_name = [SBJ_id '_latency_' stat_id '_' an_id '_' ch_list{1}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
    
    plot_data = squeeze(pk_times(:,:,1,2)); % 2nd negative peak
    violins = violinplot(squeeze(plot_data(:,:))', cond_lab, ...
        'ShowMean', true, 'ViolinAlpha', 0.3);
    
    for cond_ix = 1:numel(cond_lab)
        % Fix Mean from population estimate to sample mean
        violins(cond_ix).MeanPlot.YData = repmat(nanmean(plot_data(cond_ix,:),2), [1 2]);
        [~,mean_ix] = min(abs(violins(cond_ix).ViolinPlot.YData - nanmean(plot_data(cond_ix,:),2)));
        mean_width  = violins(cond_ix).ViolinPlot.XData(mean_ix)-cond_ix;
        violins(cond_ix).MeanPlot.XData = [cond_ix-mean_width cond_ix+mean_width];
        violins(cond_ix).MeanPlot.LineWidth = 3;
        
        % Change the colors to match condition
        %   Violin
        violins(cond_ix).ViolinColor = cond_colors{cond_ix};
        %   Box plot
        violins(cond_ix).BoxPlot.FaceColor = cond_colors{cond_ix};
        violins(cond_ix).EdgeColor = cond_colors{cond_ix};
        % 	Scatter
        %         scat_colors = zeros([numel(violins(cond_ix).ScatterPlot.XData) 3]);
        %         for s = 1:numel(SBJs)
        %             scat_colors(s,:) = cond_colors{cond_ix};
        %         end
        %         violins(cond_ix).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
        %         violins(cond_ix).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
        violins(cond_ix).ScatterPlot.CData = cond_colors{cond_ix};%scat_colors;
    end
    
    % Flip Violins horizontal
    view([90 -90]);
    
    % Add label and min RT for perspective
    ax = gca;
    ax.YLim = st.pk_lim(2,:);
    ax.YLabel.String = 'FRN Peak Latency (s)';
    ax.Title.String = ch_list{1};
    set(ax,'FontSize',16');
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        % Ensure vector graphics if saving
        if any(strcmp(fig_ftype,{'svg','eps'}))
            set(gcf, 'Renderer', 'painters');
        end
        saveas(gcf,fig_fname);
    end
end

end
