function SBJ04c_ERP_grp_stats_ANOVA(SBJs,proc_id,an_id,stat_id,save_fig,varargin)
% Compute grand average peak-to-peak FRN via jackknife:
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   conditions [str] - label of conditions to compute ERPs for
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
% OUTPUTS:
%   grp_erp [cell] - cell array with outputs of ft_timelockgrandaverage for each condition
%   NOT stat [ft struct] - output of ft_timelockstatistics, not done yet

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; app_dir = 'Users/aasthashah/Applications/fieldtrip';
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
%         elseif strcmp(varargin{v},'write_report')
%             write_report = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' stat_id '/' an_id '/'];
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

% Select Conditions of Interest
[grp_lab, ~, ~] = fn_group_label_styles(st.model_lab);
[cond_lab, cond_colors, ~, ~] = fn_condition_label_styles(st.model_lab);
% if ~strcmp(st.model_lab,{'DifOut','Out'}); error('not ready for surprise trials!'); end
grp_cond_lab = cell(size(grp_lab));
for grp_ix = 1:numel(grp_lab)
    [grp_cond_lab{grp_ix}, ~, ~, ~] = fn_condition_label_styles(grp_lab{grp_ix});
end

% Load example data
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m'];
eval(SBJ_vars_cmd);
tmp = load([root_dir 'PRJ_Error_eeg/data/',SBJs{1},'/04_proc/',SBJs{1},'_',an_id,'.mat'],'roi');
ch_list  = tmp.roi.label;
cfg = []; cfg.latency = st.stat_lim;
st_roi = ft_selectdata(cfg, tmp.roi);
time_vec = st_roi.time{1};

% Load all data
erps = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
for s = 1:length(SBJs)
    fprintf('-------------------- Processing %s ------------------------\n',SBJs{s});
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    load([SBJ_vars.dirs.proc,SBJs{s},'_',an_id,'.mat']);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    
    % Average data in stat window
    st_roi = ft_selectdata(cfg, roi);
    
    % Average data across trials
    cond_idx = fn_condition_index(cond_lab, bhv);
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for trl_ix = 1:numel(cond_trial_ix)
            trials(:,trl_ix,:) = st_roi.trial{cond_trial_ix(trl_ix)};
        end
        erps(cond_ix,s,:,:) = mean(trials,2);
    end
    
    clear tmp SBJ_vars bhv roi trial_means st_roi cond_trial_ix cond_idx
end

%% Create Design Matrix
% Design for trials + SBJ (no electrode factor yet)
design = cell(size(grp_lab));%nan([numel(cond_lab)*numel(SBJs) numel(grp_lab)]);%+1 (SBJ factor...)
for grp_ix = 1:numel(grp_lab)
    grp_assign = zeros([numel(cond_lab) 1]);
    for cond_ix = 1:numel(grp_cond_lab{grp_ix})
        grp_assign(~cellfun(@isempty,strfind(cond_lab,grp_cond_lab{grp_ix}{cond_ix}))) = cond_ix;
    end
    design{grp_ix} = repmat(grp_assign,numel(SBJs),1);
    if any(isnan(design{grp_ix})); error('NaN in ANOVA data!'); end
end
% design(:,end) = repmat([1:numel(SBJs)]',numel(cond_lab),1);

%% Prepare Data for Stats
% Jackknife procedure
if strcmp(st.grp_method,'jackknife')
    gavg = nan(size(erps));
    % Average data across sub-samples
    for s = 1:numel(SBJs)
        sbj_idx = setdiff(1:numel(SBJs),s);
        gavg(:,s,:,:) = mean(erps(:,sbj_idx,:,:),2);
    end
else
    gavg = erps;
end

% Compute data for statistics
if strcmp(st.measure,'mean')
    data = mean(gavg,4);
elseif strcmp(st.measure,'p2p')
    pk_amp = nan([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
    pk_lat = nan([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
    pk_rng = zeros([2 2]);
    for pk = 1:2
        [~, pk_rng(1,pk)] = min(abs(time_vec-st.pos_pk_lim(pk)));
        [~, pk_rng(2,pk)] = min(abs(time_vec-st.neg_pk_lim(pk)));
    end
    data = nan([numel(cond_lab) numel(SBJs) numel(ch_list)]);
    for cond_ix = 1:numel(cond_lab)
        for ch_ix = 1:numel(ch_list)
            for s = 1:numel(SBJs)
                [pk_amp(cond_ix,s,ch_ix,1), pk_lat(cond_ix,s,ch_ix,1)] = ...
                    findpeaks(squeeze(gavg(cond_ix,s,ch_ix,pk_rng(1,1):pk_rng(1,2))));
                [pk_amp(cond_ix,s,ch_ix,2), pk_lat(cond_ix,s,ch_ix,2)] = ...
                    findpeaks(-1*squeeze(gavg(cond_ix,s,ch_ix,pk_rng(2,1):pk_rng(2,2))));
                pk_amp(cond_ix,s,ch_ix,2) = -1*pk_amp(cond_ix,s,ch_ix,2);
                % Compute peak-to-peak difference
                data(cond_ix,s,ch_ix) = pk_amp(cond_ix,s,ch_ix,1) - pk_amp(cond_ix,s,ch_ix,2);
            end
        end
    end
    pk_times = zeros(size(pk_lat));
    pk_times(:,:,:,1) = time_vec(pk_lat(:,:,:,1)+pk_rng(1,1));
    pk_times(:,:,:,2) = time_vec(pk_lat(:,:,:,2)+pk_rng(2,1));
else
    error(['unknown st method for ' stat_id]);
end

% Reshape means for ANOVA
st_data = reshape(data,[numel(cond_lab)*numel(SBJs), numel(ch_list)]);
if any(isnan(st_data(:))); error('NaN in ANOVA data!'); end

%% Plot Data
if st.plot_erps
    % Plot original ERPs
    ch_figs = gobjects(size(ch_list));
    for ch_ix = 1:numel(ch_list)
        tmp = erps(:,:,ch_ix,:);
        ylims = [min(tmp(:)) max(tmp(:))];
        fig_name = ['ERPs_SBJ_' ch_list{ch_ix}];
        ch_figs(ch_ix) = figure('Name',fig_name);
        for cond_ix = 1:numel(cond_lab)
            subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix); hold on;
            title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
            plot(time_vec,squeeze(erps(cond_ix,:,ch_ix,:)),...
                'Color',cond_colors{cond_ix});%,'LineStyle',cond_styles{cond_ix});
            ylim(ylims);
        end
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
    end
    % Plot ERPs after subsampling
    if strcmp(st.grp_method,'jackknife')
        for ch_ix = 1:numel(ch_list)
            fig_name = ['ERPs_jackknife_' ch_list{ch_ix}];
            ch_figs(ch_ix) = figure('Name',fig_name);
            tmp = gavg(:,:,ch_ix,:);
            ylims = [min(tmp(:)) max(tmp(:))];
            for cond_ix = 1:numel(cond_lab)
                subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix); hold on;
                title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
                plot(time_vec,squeeze(gavg(cond_ix,:,ch_ix,:)),...
                    'Color',cond_colors{cond_ix});%,'LineStyle',cond_styles{cond_ix});
                ylim(ylims);
            end
            if save_fig
                fig_fname = [fig_dir fig_name '.' fig_ftype];
                fprintf('Saving %s\n',fig_fname);
                saveas(gcf,fig_fname);
            end
        end
    end
    % Plot peaks
    if strcmp(st.measure,'p2p')
        for ch_ix = 1:numel(ch_list)
            % ch_figs(ch_ix) = figure('Name',['ERPs_jackknife_' ch_list{ch_ix}]);
            fig_name = ['ERPs_peaks_' ch_list{ch_ix}];
            figure('Name',fig_name);
            tmp = pk_amp(:,:,ch_ix,:);
            ylims = [min(tmp(:)) max(tmp(:))];
            tmp = pk_times(:,:,ch_ix,:);
            xlims = [min(tmp(:)) max(tmp(:))];
            for cond_ix = 1:numel(cond_lab)
                subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix); hold on;
                title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
                scatter(pk_times(cond_ix,:,ch_ix,1), pk_amp(cond_ix,:,ch_ix,1), 'o', 'b');
                scatter(pk_times(cond_ix,:,ch_ix,2), pk_amp(cond_ix,:,ch_ix,2), 'o', 'r');
                for s = 1:numel(SBJs)
                    line([pk_times(cond_ix,s,ch_ix,1) pk_times(cond_ix,s,ch_ix,2)],...
                        [pk_amp(cond_ix,s,ch_ix,1) pk_amp(cond_ix,s,ch_ix,2)],...
                        'Color','k');
                end
                ylim(ylims);
                xlim(xlims);
            end
            if save_fig
                fig_fname = [fig_dir fig_name '.' fig_ftype];
                fprintf('Saving %s\n',fig_fname);
                saveas(gcf,fig_fname);
            end
        end
    end
end

%% Compute Statistics
% Omega squared isn't appropriate for interaction terms
% w2   = nan([numel(grp_lab)+1 numel(ch_list)]);
% pe2  = nan([numel(grp_lab)+1 numel(ch_list)]);
pval = nan([numel(grp_lab)+1 numel(ch_list)]);
for ch_ix = 1:numel(ch_list)
    [tmp_p, table] = anovan(st_data(:,ch_ix), design, ...
        'model', 'interaction',...% 'sstype', 2, ...% 'continuous', strmatch('RT',w2.cond),
        'varnames', grp_lab, 'display', 'off');
    
    % Adjust p for jackknife
    if strcmp(st.grp_method,'jackknife')
        f_col = strcmp(table(1,:),'F');
        df_col = strcmp(table(1,:),'d.f.');
        for grp_ix = 1:numel(grp_lab)
            grp_row = strcmp(grp_lab{grp_ix},table(:,1));
            % Adjusted F = F / (n-1)^2;
            pval(grp_ix,ch_ix) = 1-fcdf(table{grp_row,f_col}/(numel(SBJs)-1)^2, table{grp_row,df_col}, numel(SBJs)-1);
        end
        intr_row = strcmp([grp_lab{1} '*' grp_lab{2}],table(:,1));
        pval(numel(grp_lab)+1,ch_ix) = 1-fcdf(table{intr_row,f_col}/(numel(SBJs)-1)^2, table{intr_row,df_col}, numel(SBJs)-1);
    else
        pval(:,ch_ix) = tmp_p;
    end
    
    % Calculate w2 (debiased effect size; multiply with 100 to get PEV)
    %   table(:,2) = sum of squares
    %   table(:,3) = dof
    %   table(:,5) = mean square error
    %   table(:,6) = F statistic for main effects, only for rows = groups*2 (2nd set are interactions)
    %   table(:,7) = p value for main effects, only for rows = groups*2 (2nd set are interactions)
    %   rows: 1 = labels, 2:n_cond+2 = groups, n_cond+3 = DifOut interaction, end-1 = error, end = total
    %       interactions are ordered: 1*2, 1*3, ..., 1*n, 2*3, ..., 2*n, 3*4, ..., etc.
    %   w2 = (SSbw - df*MSE) / (SStot + MSE)
    %   pe2 = SSeffect / (SSeffect + SSerror)
%     mse = table{numel(grp_lab)+2,5};
%     for grp_ix = 1:numel(grp_lab)
%         grp_row = strcmp(grp_lab{grp_ix},table(:,1));
% %         w2(grp_ix,ch_ix) = (table{grp_row,2} - (table{grp_row,3} * mse))/...
% %             (table{end,2} + mse);
%     end
%     intr_row = strcmp([grp_lab{1} '*' grp_lab{2}],table(:,1));
%     w2(end,ch_ix) = (table{intr_row,2} - (table{intr_row,3} * mse))/...
%             (table{end,2} + mse);
%     clear p table
end

%% Plot results
for ch_ix = 1:numel(ch_list)
    fig_labels = {'Likely Win', 'Unlikely Loss', 'Unlikely Win', 'Likely Loss'};
    % Create and format the plot
    fig_name = ['GRP_violins_' stat_id '_' an_id '_' ch_list{ch_ix}];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
    
    violins = violinplot(squeeze(data(:,:,ch_ix))', fig_labels, ...
        'ShowMean', true, 'ViolinAlpha', 0.3);
    
    for cond_ix = 1:numel(cond_lab)
        % Fix Mean from population estimate to sample mean
        violins(cond_ix).MeanPlot.YData = repmat(mean(data(cond_ix,:,ch_ix),2), [1 2]);
        [~,mean_ix] = min(abs(violins(cond_ix).ViolinPlot.YData - mean(data(cond_ix,:,ch_ix),2)));
        mean_len = violins(cond_ix).ViolinPlot.XData(mean_ix)-cond_ix;
        violins(cond_ix).MeanPlot.XData = [cond_ix-mean_len cond_ix+mean_len];
        violins(cond_ix).MeanPlot.LineWidth = 3;
        
        % Change the colors to match condition
        %   Violin
        violins(cond_ix).ViolinColor = cond_colors{cond_ix};
        %   Box plot
        violins(cond_ix).BoxPlot.FaceColor = cond_colors{cond_ix};
        violins(cond_ix).EdgeColor = cond_colors{cond_ix};
        % 	Scatter
        scat_colors = zeros([numel(violins(cond_ix).ScatterPlot.XData) 3]);
        for s = 1:numel(SBJs)
            scat_colors(s,:) = cond_colors{cond_ix};
        end
        violins(cond_ix).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
        violins(cond_ix).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
        violins(cond_ix).ScatterPlot.CData = scat_colors;
    end
    
    % Add label and min RT for perspective
    ax = gca;
    if strcmp(st.measure,'mean')
        ax.YLabel.String = ['Mean ERP (' num2str(st.stat_lim(1)) '-' num2str(st.stat_lim(2)) ')'];
    elseif strcmp(st.measure,'p2p')
        ax.YLabel.String = 'Peak-to-Peak Amplitude';
    end
    set(ax,'FontSize',16');
    title_str = [ch_list{ch_ix} ': '];
    for grp_ix = 1:numel(grp_lab)
        title_str = [title_str grp_lab{grp_ix} ' (' num2str(pval(grp_ix,ch_ix),'%.4f') '), '];
    end
    title_str = [title_str 'Dif*Out (' num2str(pval(end,ch_ix),'%.4f') ')'];
    ax.Title.String = title_str;
    
    % Save figure
    if save_fig
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end
