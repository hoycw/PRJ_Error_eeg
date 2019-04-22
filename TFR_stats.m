function TFR_stats(SBJ, pipeline_id, an_id)

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/SCS22/','dir'); root_dir='/Users/SCS22/Desktop/Knight_Lab/';app_dir='/Users/SCS22/Documents/MATLAB/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults
%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
load([SBJ_vars.dirs.preproc SBJ '_clean_' pipeline_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' pipeline_id '_clean.mat']);
% Select Conditions of Interest
cond_easy = find(strcmp('easy', bhv.Condition));
cond_hard = find(strcmp('hard', bhv.Condition));
% Select Channel(s)
cfgs = [];
cfgs.channel = ROI;
roi = ft_selectdata(cfgs, clean_trials);
%clims = NaN([1 2]);
%comp_data = NaN([numel(roi.trial) numel(roi.time{1})]);
%for t_ix = 1:numel(roi.trial)
    %A = mean(roi.trial{t_ix},1);
    %comp_data(t_ix, :) = A;
%end
%clims(1,1) = prctile(comp_data(:),1);
%clims(1,2) = prctile(comp_data(:),99);
%clims = [min(clims(:,1)) max(clims(:,2))];
% if strcmp(cfg_tfr.method,'mtmconvol')
%     trial_lim_s_pad = [trial_lim_s(1)-max(cfg_tfr.t_ftimwin)/2 trial_lim_s(2)+max(cfg_tfr.t_ftimwin)/2+0.01];
% elseif strcmp(cfg_tfr.method,'wavelet')
%     trial_lim_s_pad = [trial_lim_s(1)-cfg_tfr.width/min(foi_center)*2 trial_lim_s(2)+cfg_tfr.width/min(foi_center)*2+0.01];
% end
%% Compute TFRs
% all cfg_tfr options are specified in the an_vars
tfr      = {};
cfg_tfr.trials = cond_easy';
tfr{1}   = ft_freqanalysis(cfg_tfr, roi);
cfg_tfr.trials = cond_hard';
tfr{2}   = ft_freqanalysis(cfg_tfr, roi);
tfr_diff = tfr{1};
tfr_diff.powspctrm = tfr{1}.powspctrm - tfr{2}.powspctrm;
%% Compute TFRs
% all cfg_tfr options are specified in the an_vars
% tfr      = {};
% for cond_ix = 1:numel(cond_lab);
%     fprintf('===================================================\n');
%     fprintf('------------- TFR Calculations for %s ----------\n',cond_lab{cond_ix});
%     fprintf('===================================================\n');
%     cfg_tfr.trials = find(cond_idx==cond_ix);
%     tfr{cond_ix}   = ft_freqanalysis(cfg_tfr, roi);
% end
%% Baseline Corrections
% do i need to do this for the diff plot
% for cond_ix = 1:2
%     switch bsln_type
%         case {'zscore', 'demean', 'my_relchange'}
%             tfr{cond_ix} = fn_bsln_ft_tfr(tfr{cond_ix},bsln_lim,bsln_type,n_boots);
%         case 'relchange'
%             cfgbsln = [];
%             cfgbsln.baseline     = bsln_lim;
%             cfgbsln.baselinetype = bsln_type;
%             cfgbsln.parameter    = 'powspctrm';
%             tfr{cond_ix} = ft_freqbaseline(cfgbsln,tfr{cond_ix});
%         otherwise
%             error(['No baseline implemented for bsln_type: ' bsln_type]);
%     end
% end

%% Stats
% Create design matrix
n_trials = numel(clean_trials.trial);
design = zeros(2,n_trials);
for an_ix = 1:2
    if an_ix==1
        design(1,(1:numel(cond_easy))) = an_ix;                                % Conditions (Independent Variable)
        design(2,(1:numel(cond_easy))) = 1:numel(cond_easy);                    % Trial Numbers
    else
        design(1,(numel(cond_easy)+1:n_trials))= an_ix; % Conditions (Independent Variable)
        design(2,(numel(cond_easy)+1:n_trials))= 1:numel(cond_hard);
    end
end

% Calculate statistics
cfg_stat.design           = design;
tfr{1}.trial = [1:numel(cond_easy)];
tfr{2}.trial = [1:numel(cond_hard)];
[stat] = ft_freqstatistics(cfg_stat, tfr{:});

%% Plot
fig_name = [SBJ 'TFA Plot'];
figure('units','normalized','Name', fig_name, 'Visible', 'on');
for cond_ix = 1:2;
    subplot(1,3, cond_ix)
    cfg = [];
    cfg.baseline     = [-0.25 0];
    cfg.baselinetype = 'absolute';
    cfg.showlabels   = 'yes';
    cfg.layout       = 'biosemi64.lay';
    cfg.ylim         = [4 20];
    %cfg.zlim = clims;
    ft_singleplotTFR(cfg, tfr{cond_ix});
    if cond_ix == 1;
        title('Easy Condition');
    else
        title('Hard Condition');
    end    
end
%% Diff Plot
    tfr_diff.statmask = stat.mask;
    cfg = [];
    cfg.baseline     = [-0.25 0];
    cfg.baselinetype = 'absolute';
    cfg.showlabels   = 'yes';
    cfg.layout       = 'biosemi64.lay';
    cfg.ylim         = [4 30];
    %cfg.zlim = clims;
    subplot(1,3, 3)
    ft_singleplotTFR(cfg, tfr_diff);
    title('Difference (Easy-Hard)')
    TFA_stack_name = [SBJ_vars.dirs.proc SBJ '_TFAPlot' '.png'];
    saveas(gcf,TFA_stack_name); 
    
    sig_ch = {};
    if sum(squeeze(stat.mask(:,:,:)))>0
        sig_ch = {sig_ch{:} stat.label{'ROI'}};
        % Link this figure to sig dir
        sig_dir = [fig_dir 'sig_ch/'];
        if ~exist(sig_dir,'dir')
            mkdir(sig_dir);
        end
        cd(sig_dir);
        link_cmd = ['ln -s ../' fig_name ' .'];
        system(link_cmd);
    end
% Save out list of channels with significant differences
sig_report_filename = [fig_dir 'sig_ch_list.txt'];
sig_report = fopen(sig_report_filename,'a');
fprintf(sig_report,'%s\n',an_id,an_id);
fprintf(sig_report,'%s\n',sig_ch{:});
fclose(sig_report);
