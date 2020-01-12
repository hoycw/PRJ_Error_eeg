function SBJ06a_CPA2(SBJ, proc_id, stat_id, plt_id)
%% Candidate-Prototype Analysis (CPA): Prototype Selection
%   Selects top IC based on spatial (elec_list) and temporal (time_win)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load processing variables
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);

%% Load the data
%loaded data from after SBJ02a --> already cleaned and trial segmented
load([SBJ_vars.dirs.preproc SBJ '_preproc_eeg_full_ft.mat']);
load([SBJ_vars.dirs.preproc SBJ '_' proc_id '_02a.mat']); %chose 02a - ica before rejection!
load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

clean_ica = ica;

[cond_lab, cond_colors, cond_styles, ~] = fn_condition_label_styles('Odd'); % maybe change this so not hardcoded
cond_idx = fn_condition_index(cond_lab, bhv);

% Create contrast: (Unexpected - Expected) for each outcome
[diff_lab, diff_pairs, diff_colors, diff_styles] = fn_condition_diff_label_styles('TarStd');

%% trim data to plotting time -- is this required?
cfg = [];
cfg.latency = cpa.time_win;
data = ft_selectdata(cfg,clean_ica);
time_vec = data.time{1};

% Get trials for plotting
trials = cell(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    cond_trial_ix = find(cond_idx==cond_ix);
    trials{cond_ix} = nan([numel(data.label) numel(cond_trial_ix) numel(data.time{1})]);
    for t_ix = 1:numel(cond_trial_ix)
        trials{cond_ix}(:,t_ix,:) = data.trial{cond_trial_ix(t_ix)};
    end
end

%% Compute plotting data
diff_waves = zeros(numel(data.label), numel(diff_lab), numel(data.time{1}));
[~, min_ix] = min(abs(clean_ica.time{1,1}(:) - time_win(1)));
[~, max_ix] = min(abs(clean_ica.time{1,1}(:) - time_win(2)));
for comp_ix = 1:numel(data.label)
    % Compute means and variance
    means = NaN([numel(cond_lab) numel(data.time{1})]);
    sems  = NaN([numel(cond_lab) numel(data.time{1})]);
    for cond_ix = 1:numel(cond_lab)
        means(cond_ix,:) = squeeze(mean(trials{cond_ix}(comp_ix,:,:),2));
        sems(cond_ix,:) = squeeze(std(trials{cond_ix}(comp_ix,:,:),[],2))./sqrt(size(trials{cond_ix},2))';
        means_all(comp_ix, cond_ix, :) = means(cond_ix,:);
        sems_all(comp_ix, cond_ix, :) = sems(cond_ix,:);
    end
    
    plot_means = NaN([numel(diff_lab) numel(data.time{1})]);
    % Compute stats between conditions in contrast within time windows
    %!!! initialize variables
    for diff_ix = 1:numel(diff_lab)
        %!!! don't loop over all time, only cpa_lim
        for time_ix = 1:numel(data.time{1})
            [sig(comp_ix, time_ix), p(comp_ix, time_ix)] = ttest2(trials{diff_pairs{diff_ix}(1)}(comp_ix, :, time_ix), trials{diff_pairs{diff_ix}(2)}(comp_ix, :, time_ix));
        end
        sig_time_win(comp_ix,:) = sig(comp_ix, min_ix: max_ix);
        total_num_sig(comp_ix) = length(find(sig_time_win(comp_ix,:) == 1));
        fraction_sig(comp_ix) = total_num_sig(comp_ix)/(max_ix - min_ix + 1);
        lenmax = 1;
        len = 1;
        %find longest stretch of significance
        for n = 2:numel(sig_time_win(comp_ix,:))
            if sig_time_win(comp_ix, n) == sig_time_win(comp_ix, n-1) && sig_time_win(comp_ix, n) == 1
                len = len+1;
            else
                if len > lenmax
                    lenmax = len;
                end
                len = 1;
            end
        end
        if len > lenmax
                lenmax = len;
        end
        sig_length_max(comp_ix) = lenmax;
        %Difference Wave (for debugging)
        plot_means(diff_ix,:) = means(diff_pairs{diff_ix}(1),:)-means(diff_pairs{diff_ix}(2),:);
    end
    diff_waves(comp_ix, diff_ix,:) = plot_means(diff_ix,:); % whole time period
end

% Select components with consecutive significance
[~, erp_components] = find(sig_length_max > cpa.min_sig_len);

%% compute ERP for each Electrode
cfg = [];
cfg.channel = 'all';
erps = ft_timelockanalysis(cfg, data);
num_elec = floor(numel(clean_ica.topo(1,:))/12);
max_elecs = zeros(numel(clean_ica.topo(1,:)), num_elec);
%Find max elecs for each component
%!!! initialize
for comp_ix = 1:numel(clean_ica.topo(1,:))
    [~, max_elecs(comp_ix,:)] = maxk(abs(clean_ica.topo(:,comp_ix)), num_elec);
    for max_elecs_ix = 1:num_elec
        max_names(comp_ix,max_elecs_ix) = clean_ica.topolabel(max_elecs(comp_ix,max_elecs_ix));
    end
end

% Find the index of the target electrodes
for elec_ix = 1:numel(elec_list)
    temp_elec = elec_list{elec_ix};
    for topo_elec_ix = 1:numel(clean_ica.topolabel)
        temp_ica = clean_ica.topolabel{topo_elec_ix};
        if (strcmp(temp_elec, temp_ica))
            index_electrodes(elec_ix) = topo_elec_ix;
        end
    end
end
topo_components = [];
% Find components which have 2+ of the target electrodes in their max electrodes
temp_array = zeros(1, numel(data.label));
for elec_ix = 1:numel(index_electrodes)
    [temp,~] = find(max_elecs == index_electrodes(elec_ix));
    temp_array(temp) = temp_array(temp) + 1;
end
topo_components = find(temp_array > 1); 
topo_components = topo_components';
components = intersect(topo_components, erp_components);

% Find ICs with most relative explained variance
%   WARNING: I did this by relative explanation of variance. Does not use significance;
%   can't find a way to find the relative explanation of variance with the ICA data.
%   Works okay with out this but pulls extra components.
components_final_ix = find(components < cpa.ic_rank_max);
components_final = components(components_final_ix);
disp(components_final);

%% PLOT AND SAVE
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/CPA_Graphs/' SBJ '/plot/'] ;
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

sig_chunks = cell(numel(data.label), 1);
for comp_ix = 1: numel(data.label)
    target_elec_lab = strjoin(cpa.elec_list,'_');
    fig_name = [SBJ 'Component0' num2str(comp_ix) 'Electrodes_' target_elec_lab '_Oddball_'];
    figure('Name', fig_name, 'units','normalized',...
        'outerposition',[0 0 0.5 0.8], 'Visible', 'off');
    ebars = cell(size(cond_lab));
    main_lines = gobjects([numel(cond_lab)+1 1]);
    for cond_ix = 1:numel(cond_lab)
        ebars{cond_ix} = shadedErrorBar(time_vec(:), means_all(comp_ix, cond_ix, :), sems_all(comp_ix, cond_ix, :),...
            {'Color',cond_colors{cond_ix},'LineWidth',plt.mean_width,...
            'LineStyle',cond_styles{cond_ix}},plt.errbar_alpha);
        hold on
        main_lines(cond_ix) = ebars{cond_ix}.mainLine;
    end
    sig_chunks{comp_ix} = fn_find_chunks(squeeze(sig_time_win(comp_ix,:)));
    sig_chunks{comp_ix}(squeeze(sig_time_win(comp_ix,sig_chunks{comp_ix}(:,1))==0),:) = [];
    data_lim = [min(min(means_all(comp_ix,:,:)-sems_all(comp_ix,:,:))) max(max(means_all(comp_ix,:,:)+sems_all(comp_ix,:,:)))];
    % Plot Significance
    for sig_ix = 1:size(sig_chunks{comp_ix},1)
        time_vec = clean_ica.time{1,1};
        sig_times = time_vec([sig_chunks{comp_ix}(sig_ix,1):sig_chunks{comp_ix}(sig_ix,2)]+min_ix-1);
        sig_y = data_lim(1) + cond_ix*data_lim(1)*plt.sig_loc_factor;
        sig_line = line(sig_times,repmat(sig_y,size(sig_times)),...
            'LineWidth',plt.sig_width,'Color',cond_colors{cond_ix});
        if sig_ix==1
            main_lines(end+1) = sig_line;
        end
    end
    % Axes and Labels
    ax.YLabel.String = 'uV';
    ax.XLim          = [plt.plt_lim(1) plt.plt_lim(2)];
    ax.XTick         = plt.plt_lim(1):plt.x_step_sz:plt.plt_lim(2);
    ax.XLabel.String = 'Time (s)';
    target_elec_lab = '';
    for elec_ix = 1:3
        target_elec_lab = [target_elec_lab max_names{comp_ix,elec_ix} ' '];
    end
    title(['Electrodes: ' target_elec_lab ' Time Window: ' num2str((sig_length_max(comp_ix)-1)*4) 'ms Percent Significant: ' num2str(fraction_sig(comp_ix)*100)]); %4 is used because thats the number of miliseconds in 1 sample
    leg_lab = [cond_lab 'F' cond_lab(cond_ix)];
    %if plt.legend
       legend([main_lines{:}],leg_lab,'Location',plt.legend_loc);
   % end
     fig_fname = [fig_dir fig_name '.png'];
     fprintf('Saving %s\n',fig_fname);
     % Ensure vector graphics if saving
     %saveas(gcf,fig_fname);
     
     %This copies the files over to the correct folder so it can be
     %compared (these figures are generated in 02a)
     comp_label = clean_ica.label(comp_ix,1);
     comp_label = comp_label{1};
     fig_dir_odd = [SBJ_vars.dirs.proc_stack SBJ comp_label '_ICA_plots_odd.png'];
     stack_dir_odd = [SBJ_vars.dirs.proc_stack SBJ comp_label '_ERP_stack_odd.png'];
     fig_dir_tt = [SBJ_vars.dirs.proc_stack SBJ comp_label '_ICA_plots.png'];
     stack_dir_tt = [SBJ_vars.dirs.proc_stack SBJ comp_label '_ERP_stack.png'];
     if exist(fig_dir_odd, 'file')
         copyfile(fig_dir_odd, fig_dir);
     end
     if exist(stack_dir_odd, 'file')
         copyfile(stack_dir_odd, fig_dir);
     end
     if exist(stack_dir_tt, 'file')
         copyfile(stack_dir_tt, fig_dir);
     end
     if exist(fig_dir_tt, 'file')
         copyfile(fig_dir_tt, fig_dir);
     end
    
end   


%% Save Data
clean_data_fname = [SBJ_vars.dirs.preproc SBJ '_' proc_id '_06a.mat'];
save(clean_data_fname, '-v7.3', 'components', 'components_final');

end
