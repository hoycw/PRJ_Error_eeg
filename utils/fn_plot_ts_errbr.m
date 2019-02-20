function fn_plot_ts_errbr(plot_info, data, var, event, cond)
% Plot mean time series with shaded error bars on given plot handle
% INPUTS:
%   plot_info [struct] - contains general plotting details
%       .fig [figure handle] - handle of the figure where the plot should appear
%       .ax [axis handle] - handle of the subplot/axis of fig where the plot should appear
%       .title [str] - title for the plot
%       .x_step [int] - step size for x axis tick marks
%       .x_lab [int array] - list of labels for x axis tick marks
%       .legend [0/1] - flag to plot legend
%       .legend_loc [str] - location of legend
%   data [matrix] - mean traces to be plotted; rows = time series already trimmed to plot limits
%   var [matrix] - variance around mean traces for plotting error bars (also pre-trimmed)
%   event [struct] - contains details of any events to be plotted
%       .time [double] - index of events for line plot (e.g., stimulus or response)
%       .name [str] - name for legend
%       .width [double] - width of line
%       .color [str] - color to plot event
%       .style [str] - line style to plot
%   cond [struct] - contains details of the conditions
%       .name [cell str array] - strings that are names of the time series for legend
%                                   leave empty if no legend desired
%       .style [cell str array] - strings indicating line styles for each time series
%       .color [cell str array] - strings of colors for each time series
%       .alpha [double] - transparency of error bars
% TO BE ADDED:
%   sig [matrix] - list of time point with significant differences


% Set Up Parameteres
set(0, 'currentfigure', plot_info.fig);   % set fig as current figure
set(plot_info.fig, 'currentaxes', plot_info.ax);    % set ax as current axes in figure fig

% Plot Event Locked Average
hold on;
ebars = {};
main_lines = [];
for cond_ix = 1:length(cond.name)
    ebars{cond_ix} = shadedErrorBar([], data(cond_ix,:), var(cond_ix,:),...
        {'Color',[cond.color{cond_ix}],...
        'LineStyle',cond.style{cond_ix}},cond.alpha(cond_ix));
    main_lines = [main_lines ebars{cond_ix}.mainLine];
end
% ylim(ylims);
ylims = ylim;

% Plot Extra Features (events, significance)
for evnt_ix = 1:numel(event.name)
    event_line = line([event.time(evnt_ix) event.time(evnt_ix)],ylims,...
        'LineWidth',event.width(evnt_ix),'Color',event.color{evnt_ix},'LineStyle',event.style{evnt_ix});
    main_lines = [main_lines event_line];
end
% % Significance marker
% scatter(sig_times{an_ix}{ch_ix}-plot_lim(1),zeros([1 length(sig_times{an_ix}{ch_ix})]),...
%     'Marker','.','LineWidth',0.1,'MarkerEdgeColor','k');


% Axes and Labels
plot_info.ax.YLim          = ylims;
plot_info.ax.YLabel.String = 'HFA (z-score)';
plot_info.ax.XLim          = [0,size(data,2)];
plot_info.ax.XTick         = 0:plot_info.x_step:size(data,2);
plot_info.ax.XTickLabel    = plot_info.x_lab;
plot_info.ax.XLabel.String = 'Time (s)';
title(plot_info.title);
if plot_info.legend==1
    legend(main_lines,cond.name{:},event.name{:},'Location',plot_info.legend_loc);
end

%         % Plot Single Trials Per Condition
%         for cond_ix = 1:length(cond_lab)
%             subplot(length(analyses),length(cond_lab)+1,...
%                 an_ix*(length(cond_lab)+1)-length(cond_lab)+cond_ix); hold on;
%
%             imagesc(squeeze(trials(an_ix,ch_ix,RTs_sort_idx(cond_idx(cond_ix,:)),win_on:win_off)));
%             set(gca,'YDir','normal');
%             if strcmp(event,'stim')
%                 scat = scatter(RTs_sort_cond{cond_ix}+plot_lim(1),1:sum(cond_idx(cond_ix,:)),...
%                     'MarkerFaceColor',[cond_colors{cond_ix}],...
%                     'MarkerEdgeColor','k');
%             end
%             ylim([1 sum(cond_idx(cond_ix,:))]);
%             event_line = line([plot_lim(1) plot_lim(1)],ylim,...
%                 'LineWidth',event_ln_width,'Color','k');
%             % Axes and Labels
%             cbar = colorbar;
%             caxis(clims);
%             ax = gca;
%             ax.XLim = [0,win_off-win_on];
%             ax.XTick = 0:x_step:win_off;
%             ax.XTickLabel = x_lab;
%
%             title(strcat(header_ecog.channel_labels(ch_ix), ':', analysis_id{an_ix},...
%                 ',',cond_lab{cond_ix},' (Single Trial, ',event,'-locked)'));
%         end

end
