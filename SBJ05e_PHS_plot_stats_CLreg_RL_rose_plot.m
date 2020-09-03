function SBJ05e_PHS_plot_stats_CLreg_RL_rose_plot(SBJ_id,proc_id,an_id,stat_id,save_fig)
error('This was only used to investigate the flip in sign of CLreg coefficients');
%% Plot rose plot of phase values at both max and min TFR tiles
%   Only for single channel right now...
%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Documents/MATLAB/';
else; root_dir='/Volumes/hoycw_clust/'; app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Handle Variable Inputs & Defaults

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);
SBJs = fn_load_SBJ_list(SBJ_id);
%load SBJ05d data -- behavior, beta weights
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'betas','model','data');
load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/04_proc/' SBJs{1} '_' proc_id '_' an_id '.mat']);

% Load time data
cfgs = []; cfgs.latency = st.stat_lim;
st_tfr = ft_selectdata(cfgs, tfr);
st_time_vec = st_tfr.time;
%% Find TFR values

matrix_tfr_tiles = zeros(2, 2); %Dimensions -- tfr tiles (max and min), indices (time_vec, freq)
[temp_val, temp] = max(squeeze(betas(:, :,:)),[], 1); % Find Max Value for Each condition
[temp_val2, temp2] = max((temp_val),[], 2); % Find Max Value for Each Frequency
[~, matrix_tfr_tiles(1, 1)] = max(temp_val2); % Find index of Max Frequency Value
matrix_tfr_tiles(1, 2) = temp2(matrix_tfr_tiles(1, 1)); % Find the time_vec indices where that max value occurs
[temp_val, temp] = min(squeeze(betas(:, :,:)),[], 1); %Repeat for min values
[temp_val2, temp2] = min((temp_val),[], 2); % Find Max Value for Each Frequency
[~, matrix_tfr_tiles(2, 1)] = min(temp_val2);
matrix_tfr_tiles(2, 2) = temp2(matrix_tfr_tiles(2, 1));
%{
for cond_ix = 1:size(betas,1)
        [temp_val, temp] = max(squeeze(betas(cond_ix, :,:)),[], 2); % Find Max Value for Each Frequency
        [~, matrix_tfr_tiles(cond_ix, 1, 1)] = max(temp_val); % Find index of Max Frequency Value
        matrix_tfr_tiles(cond_ix, 1, 2) = temp(matrix_tfr_tiles(cond_ix, 1, 1)); % Find the time_vec indices where that max value occurs
        [temp_val, temp] = min(squeeze(betas(cond_ix, :,:)),[], 2); %Repeat for min values
        [~, matrix_tfr_tiles(cond_ix, 2, 1)] = min(temp_val);
        matrix_tfr_tiles(cond_ix, 2, 2) = temp(matrix_tfr_tiles(cond_ix, 2, 1));
end
%}
matrix_angles = cell(2, 2); %Dimensions -- tfr tiles (max and min), PE type (neg or pos)
matrix_mean_vector = zeros(2, 2, 2); %Dimensions -- conditions (pWin, sPE, uPE), tfr tiles (max and min), PE type (neg or pos), (theta, rho)
PE_trials{1} = find(squeeze(model(:,2)) < 0); % For each condition, find positive prediction error trials
PE_trials{2} = find(squeeze(model(:,2)) >= 0);% For each condition, find negative prediction error trials
    for tfr_ix = 1:2
        for pe_ix = 1:2
            matrix_angles(tfr_ix, pe_ix) = {data(PE_trials{pe_ix}, matrix_tfr_tiles(tfr_ix, 2), matrix_tfr_tiles(tfr_ix, 1))};
            %Take ITPC data from given conditions
            matrix_mean_vector(tfr_ix, pe_ix, 1) = circ_mean(matrix_angles{tfr_ix, pe_ix});
            % Find mean vector length and rho value
            matrix_mean_vector(tfr_ix, pe_ix, 2) = circ_r(matrix_angles{tfr_ix, pe_ix});
        end
    end

titles = {'Negative PEs: Max Beta Tile',  'Positive PEs: Max Beta Tile', 'Negative PEs: Min Beta Tile','Positive PEs: Min Beta Tile'};
    fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' an_id '/' stat_id '/'];
    fig_name = [SBJ_id '_' stat_id '_' an_id '_Rose_Plot'];
    fig_temp = figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    for tfr_ix = 1:2
    %subplot(2,2,(2*tfr_ix+1)-2); % ADDED
    %pe_colors = {[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330]};
        for pe_ix = 1:2
            subplot(2,2,(2*tfr_ix+pe_ix)-2);
            polarhistogram(matrix_angles{tfr_ix, pe_ix},[-pi:pi/5:pi],'Normalization','probability');
            fig_temp = subplot(2,2,(2*tfr_ix+pe_ix)-2);
            rL = rlim(fig_temp);
            rho_norm = matrix_mean_vector(tfr_ix, pe_ix, 2) * range(rL) + rL(1);
            hold on
            polarplot([matrix_mean_vector(tfr_ix, pe_ix, 1) matrix_mean_vector(tfr_ix, pe_ix, 1)], [0 rho_norm], 'Color', 'k', 'LineWidth', 2);
            title([titles{(2*tfr_ix+pe_ix)-2} ' Frequency: '  num2str(matrix_tfr_tiles(tfr_ix, 2)) ' Hz, Time: ' num2str(st_time_vec(matrix_tfr_tiles(tfr_ix, 1))) ' sec']); 
            %{
            polarhistogram(matrix_angles{cond_ix, tfr_ix, pe_ix},[-pi:pi/5:pi],'Normalization','probability', 'FaceColor', pe_colors{pe_ix}, 'FaceAlpha', 0.5);
            fig_temp = subplot(2,2,(2*tfr_ix+1)-2);
            rL = rlim(fig_temp);
            rho_norm = matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 2) * range(rL) + rL(1);
            hold on
            polarplot([matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 1) matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 1)], [0 rho_norm], 'LineWidth', 2, 'Color', pe_colors{pe_ix});
            title([titles{(2*tfr_ix+pe_ix)-2} ' Frequency: '  num2str(matrix_tfr_tiles(cond_ix, tfr_ix, 1)) ' Hz, Time: ' num2str(st_time_vec(matrix_tfr_tiles(cond_ix, tfr_ix, 2))) ' sec']);
            hold on
            %}
        end
    %end
    if save_fig
        fig_filename = [fig_dir fig_name '.png'];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end
end
