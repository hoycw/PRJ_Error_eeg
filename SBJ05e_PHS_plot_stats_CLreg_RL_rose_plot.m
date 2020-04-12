function SBJ05e_PHS_plot_stats_CLreg_RL_rose_plot(SBJ_id,proc_id,an_id,stat_id,save_fig)
%% Plot rose plot of phase values at both max and min TFR tiles
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

% Define default options
if ~exist('fig_vis','var'); fig_vis = 'on'; end
if ~exist('fig_ftype','var'); fig_ftype = 'png'; end
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Load Results
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

%load SBJ05d data -- behavior, beta weights
load([root_dir 'PRJ_Error_eeg/data/GRP/' SBJ_id '_' stat_id '_' an_id '.mat'],'betas','model','data');

%% Find TFR values
%{
max_time_vec = zeros(3, numel);
i_max_time_vec = zeros(1, 3);
max_freq = zeros(1, 3);
i_max_freq = zeros(1, 3);
max_time_vec = zeros(1, 3);
for ix = 1:numel(reg_lab)
    [max_time_vec(ix), i_max_time_vec(ix)] = max(squeeze(betas(ix, :,:)),[], 2);
    [max_freq(ix), i_max_freq(ix)] = max(max_time_vec(ix));
    i_max_time_vec(ix) = i_max_time_vec(ix, i_max_freq);
    [min_time_vec(ix), i_min_time_vec(ix)] = min(squeeze(betas(ix, :,:)),[], 2);
    [min_freq(ix), i_min_freq(ix)] = min(min_time_vec(ix));
    i_min_time_vec(ix) = i_min_time_vec(ix, i_min_freq);
end  
for ix = 1:numel(reg_lab)
    neg_PE_trials(ix) = find(squeeze(model(ix,:)) < 0);
    pos_PE_trials(ix) = find(squeeze(model(ix,:)) >= 0);
    max_angles_neg(ix) = data(neg_PE_trials(ix), i_freq, i_max_time_vec);
    min_angles_neg(ix) = data(neg_PE_trials(ix),i_freq, i_min_time_vec);
    max_angles_pos(ix) = data(pos_PE_trials(ix), i_freq, i_max_time_vec);
    min_angles_pos(ix) = data(pos_PE_trials(ix),i_freq, i_min_time_vec);
end
%}
matrix_tfr_tiles = zeros(numel(betas,1), 2, 2); %Dimensions -- conditions (pWin, sPE, uPE), tfr tiles (max and min), indices (freq, time_vec)
for cond_ix = 1:size(betas,1)
        [temp_val, temp] = max(squeeze(betas(cond_ix, :,:)),[], 2); % Find Max Value for Each Frequency
        [~, matrix_tfr_tiles(cond_ix, 1, 1)] = max(temp_val); % Find index of Max Frequency Value
        matrix_tfr_tiles(cond_ix, 1, 2) = temp(matrix_tfr_tiles(cond_ix, 1, 1)); % Find the time_vec indices where that max value occurs
        [temp_val, temp] = min(squeeze(betas(cond_ix, :,:)),[], 2); %Repeat for min values
        [~, matrix_tfr_tiles(cond_ix, 2, 1)] = min(temp_val);
        matrix_tfr_tiles(cond_ix, 2, 2) = temp(matrix_tfr_tiles(cond_ix, 2, 1));
end        
matrix_angles = cell(numel(betas,1), 2, 2); %Dimensions -- conditions (pWin, sPE, uPE), tfr tiles (max and min), PE type (neg or pos)
matrix_mean_vector = zeros(numel(betas, 1), 2, 2, 2); %Dimensions -- conditions (pWin, sPE, uPE), tfr tiles (max and min), PE type (neg or pos), (theta, rho)
for cond_ix = 1:size(betas, 1)
    PE_trials{1} = find(squeeze(model(:,cond_ix)) < 0); % For each condition, find positive prediction error trials
    PE_trials{2} = find(squeeze(model(:,cond_ix)) >= 0);% For each condition, find negative prediction error trials
    for tfr_ix = 1:2
        for pe_ix = 1:2
            matrix_angles(cond_ix, tfr_ix, pe_ix) = {data(PE_trials{pe_ix}, matrix_tfr_tiles(cond_ix, tfr_ix, 1), matrix_tfr_tiles(cond_ix, tfr_ix, 2))};
            %Take ITPC data from given conditions
            matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 1) =circ_mean(matrix_angles{cond_ix, tfr_ix, pe_ix});
            % Find mean vector length and rho value
            matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 2) = circ_r(matrix_angles{cond_ix, tfr_ix, pe_ix});
        end
    end
end
cond_lab = {'PWin', 'sPE', 'uPE'};
titles = {'Negative Prediction Errors -- Max TFA Tile Phase Values', 'Negative Prediction Errors -- Min TFA Tile Phase Values'...
    'Positive Prediction Errors -- Max TFA Tile Phase Values', 'Positive Prediction Errors -- Min TFA Tile Phase Values'};
for cond_ix = 1:size(betas,1)
    fig_dir = [root_dir 'PRJ_Error_eeg/results/TFR/' an_id '/' stat_id '/'];
    fig_name = [SBJ_id '_' stat_id '_' an_id '_Rose_Plot_' cond_lab{cond_ix}];
    fig_temp = figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 0.8 0.8],'Visible',fig_vis);
    for tfr_ix = 1:2
        for pe_ix = 1:2
            subplot(2,2,(2*tfr_ix+pe_ix)-2);
            polarhistogram(matrix_angles{cond_ix, tfr_ix, pe_ix},[-pi:pi/5:pi],'Normalization','probability');
            fig_temp = subplot(2,2,(2*tfr_ix+pe_ix)-2);
            rL = rlim(fig_temp);
            rho_norm = matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 2) * range(rL) + rL(1);
            hold on
            polarplot([matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 1) matrix_mean_vector(cond_ix, tfr_ix, pe_ix, 1)], [0 rho_norm], 'Color', 'k');
            title([titles{(2*tfr_ix+pe_ix)-2} ' Frequency: '  num2str(matrix_tfr_tiles(cond_ix, tfr_ix, 1)) ' Time Vec: ' num2str(matrix_tfr_tiles(cond_ix, tfr_ix, 2))]);
        end
    end
    if save_fig
        fig_filename = [fig_dir fig_name '.png'];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
    end
end
end
