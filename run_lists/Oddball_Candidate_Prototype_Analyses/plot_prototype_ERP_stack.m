for f_ix = 1:numel(final_ics)
    comp_ix = final_ics(f_ix);
    %% Compute plotting data    
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(clean_ica.time{1})]);
    sems  = NaN([numel(cond_lab) numel(clean_ica.time{1})]);
    stack = [];%NaN([numel(cond_idx) numel(clean_ica.time{1})]);
    trl_cond = zeros(size(cond_idx)); trl_ix = 0;
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(comp_ix,:,:),2));
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(comp_ix,:,:),[],2))./sqrt(size(trials{cond_ix},2))';
        
        stack = vertcat(stack,squeeze(trials{cond_ix}(comp_ix,:,:)));
        trl_ix = trl_ix + 1;
        trl_cond(trl_ix:trl_ix+size(trials{cond_ix},2)-1) = repmat(cond_ix,[size(trials{cond_ix},2) 1]);
        trl_ix = trl_ix+size(trials{cond_ix},2)-1;
    end
    
    %% Create plot
    fig_name = [SBJ '_' cpa_id '_' num2str(final_ics(f_ix))];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 1],'Visible',fig_vis);
    
    %% Plot Trial Stack
    axes = gobjects([2 1]);
    axes(1) = subplot(3,1,[1 2]); hold on;
    
    imagesc(clean_ica.time{1},1:numel(clean_ica.trial),stack);
    set(gca,'YDir','Normal');
    
    % Plot Events
    for evnt_ix = 1:numel(plt.evnt_lab)
        if any(strcmp(plt.evnt_lab{evnt_ix},{'F','Fon'}))
            for trl_ix = 1:numel(clean_ica.trial)
                scatter(evnt_times(evnt_ix), trl_ix,...
                    'MarkerFaceColor',[cond_colors{trl_cond(trl_ix)}],'MarkerEdgeColor','none',...
                    'Marker',cond_markers{trl_cond(trl_ix)});
            end
        else
            line([evnt_times(evnt_ix) evnt_times(evnt_ix)],[1 numel(clean_ica.trial)],...
                'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
                'LineStyle',plt.evnt_style);
        end
    end
    
    axes(1).YLabel.String = 'Trials';
    axes(1).YLim = [1 numel(clean_ica.trial)];
    axes(1).XLim = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(1).XLabel.String = 'Time (s)';
    title(num2str(comp_ix));
%     if plt.legend
%         legend(an.event_type,'Location',plt.legend_loc);
%     end
    colorbar('Location','northoutside');
    
    %% Plot ERPs
    axes(2) = subplot(3,1,3); hold on;
    
    % Plot individual trials per condition
    if plt.butterfly
        for cond_ix = 1:numel(cond_lab)
            plot(clean_ica.time{1},squeeze(trials{cond_ix}(ch_ix,:,:)),...
                'Color',cond_colors{cond_ix},'LineWidth',plt.butterfly_width,...
                'LineStyle',cond_styles{cond_ix});
        end
    end
    
    % Plot Means (and variance)
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+numel(plt.evnt_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(clean_ica.time{1}, means(cond_ix,:), sems(cond_ix,:),...
            'lineProps',{'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},'patchSaturation',plt.errbar_alpha);
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    
    % Plot Extra Features (events, significance)
    for evnt_ix = 1:numel(plt.evnt_lab)
        main_lines(numel(cond_lab)+evnt_ix) = line([evnt_times(evnt_ix) evnt_times(evnt_ix)],ylim,...
            'LineWidth',plt.evnt_width,'Color',plt.evnt_color,...
            'LineStyle',plt.evnt_styles{evnt_ix});
    end
    leg_lab = [cond_lab an.event_type];
    
    % Axes and Labels
    axes(2).YLabel.String = 'uV';
    axes(2).XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    axes(2).XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    axes(2).XLabel.String = 'Time (s)';
    axes(2).Title.String  = [clean_ica.label{comp_ix} ': ' conditions];
    if plt.legend
        legend(main_lines,leg_lab{:},'Location','best');%plt.legend_loc);
    end
    
end
