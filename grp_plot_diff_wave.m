function grp_plot_diff_wave(SBJs, plt_id, an_id, fig_vis, save_fig, fig_ftype, conds)
%% Check which root directory
if exist('/home/knight/','dir');root_dir='/home/knight/';ft_dir=[root_dir 'PRJ_Error_eeg/Apps/fieldtrip/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';ft_dir='/Users/sheilasteiner/Downloads/fieldtrip/';
elseif exist ('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/fieldtrip-private']);
addpath(ft_dir);
ft_defaults

%% Load preprocessed data
for s = 1:length(SBJs)
    SBJ_vars_cmd{s} = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd{s});
    SBJ_vars_all{s} = SBJ_vars;
end
plt_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
cond_lab_hit = [1 0];
cond_lab = {'easy', 'hard'};
cond_colors = {[0.6350, 0.0780, 0.1840], [0.3010, 0.7450, 0.9330], 	[0.4660, 0.6740, 0.1880], [0.4940, 0.1840, 0.5560]};
for s = 1:length(SBJs)
    clean_stats_fname{s} = [SBJ_vars_all{s}.dirs.proc SBJs{s} '_' an_id '.mat'];
    load(clean_stats_fname{s});
end

% load all SBJ roi_erp
roi_erp = cell(size(SBJs));
for s = 1:length(SBJs)
    tmp = load(clean_stats_fname{s}); roi_erp{s} = tmp.roi_erp;
end

%% Plot ERPs and stats
if save_fig
    fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/FRN/diffeence_' an_id conds '/'];
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
end

fig_name = ['FRN_dif' conds];
f = figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 0.5 0.5],'Visible', fig_vis);
if strcmp(conds, 'DifWL');
    index_wav_1 = [1 3];
    index_wav_2 = [4 2];
    cond_info.name       = {'Easy Win - Hard Win', 'Hard Loss - Easy Loss'};
elseif strcmp(conds, 'DifEH');
    index_wav_1 = [1 2];
    index_wav_2 = [4 3];
    cond_info.name  = {'Easy Win - Easy Loss', 'Hard Loss - Hard Win'};
else
    sprintf('Unknown condition variable!');
end
for s = 1:length(SBJs)
    diff_wave = roi_erp;
    diff_wave{s}{1}.avg = roi_erp{s}{index_wav_1(1)}.avg - roi_erp{s}{index_wav_1(2)}.avg;
    diff_wave{s}{2}.avg = roi_erp{s}{index_wav_2(1)}.avg - roi_erp{s}{index_wav_2(2)}.avg;
    rmfield(diff_wave{s}{1}, 'trial');
    rmfield(diff_wave{s}{2}, 'trial');
    rmfield(diff_wave{s}{2}, 'var');
    rmfield(diff_wave{s}{1}, 'var');
    rmfield(diff_wave{s}{2}, 'var');
end

gstat = cell(size(diff_wave));
averages = cell(numel(SBJs));
cfg_avg.paramater = 'avg';
cfg_avg.keepindividual = 'no';

for cond_ix = 1:4
    for index = 1: numel(SBJs)
        averages{index} = diff_wave{s}{cond_ix};
    end
    gstat{cond_ix} = ft_timelockgrandaverage(cfg_avg, averages{1:numel(SBJs)});
end

 plot_info.fig        = f;
    plot_info.sig_color  = plt_vars.sig_color;
    plot_info.sig_alpha  = plt_vars.sig_alpha;
    plot_info.x_step = plt_vars.x_step_sz*plt_vars.fsample;
    plot_info.x_lab  = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    plot_info.y_lab  = 'uV';
    plot_info.legend = plt_vars.legend;
    plot_info.legend_loc = plt_vars.legend_loc;
    % Event plotting params
    event_info.name      = {an.event_type};
    %%
    [~,event_info.time]  = min(abs(roi_erp{1}{1}.time-1.8));
    event_info.width     = plt_vars.evnt_width;
    event_info.color     = {plt_vars.evnt_color};
    event_info.style     = {plt_vars.evnt_style};
    % Condition plotting params
    cond_info.style      = {'-', '-', '-', '-'};
    cond_info.color      = cond_colors;
    cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 4]);

    %% Plot all ERPs together
    subplot(numel(cond_lab),1,1); hold on;
    plot_info.title  = 'Difference Waves';
    plot_info.ax     = gca;
%% calculate statistics
        % Compute means and variance
      
    
        means = NaN([numel(cond_lab)+numel(cond_lab_hit) size(gstat{1}.avg,2)]);
        sem   = NaN([numel(cond_lab)+numel(cond_lab_hit) size(gstat{1}.avg,2)]);
        temp = cell(size(numel(SBJs)));
        store = zeros(numel(SBJs),length(gstat{1}.avg));
        for c_ix = 1:numel(cond_lab)+numel(cond_lab_hit)
            means(c_ix,:) = gstat{c_ix}.avg(1,:);
            for s = 1:numel(SBJs)
                temp{s} = roi_erp{s}{c_ix};
                for n = 1: length(gstat{1}.avg)
                    store(s,n) = temp{s}.avg(n);
                end
            end
            for t_ix = 1: length(gstat{1}.avg)
                sem(c_ix,t_ix) = std(store(:,t_ix))/ sqrt(numel(SBJs));
            end
        end
        
    
    fn_plot_ts_errbr(plot_info,means,sem,event_info,cond_info);
        
    