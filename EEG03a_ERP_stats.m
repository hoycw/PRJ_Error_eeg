function EEG03a_ERP_stats(SBJ,conditions,pipeline_id,an_id)
% Calculates ERPs, computes cluster-based statistics, and plots the results
% clear all; %close all;

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
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

% Load Data
load([SBJ_vars.dirs.preproc SBJ '_clean_' pipeline_id '.mat']);
load([SBJ_vars.dirs.events SBJ '_behav_' pipeline_id '_clean.mat']);

% Select Conditions of Interest
[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);
cond_idx = fn_condition_index(conditions, bhv);

% Select Channel(s)
cfgs = [];
cfgs.channel = ROI;
roi = ft_selectdata(cfgs, clean_trials);

%% Compute ERPs
% Preprocess the data
cfgpp = [];
cfgpp.hpfilter     = hp_yn;
cfgpp.hpfreq       = hp_freq;
cfgpp.hpfiltord    = hp_filtord;            % Leaving blank causes instability error, 1 or 2 works 
cfgpp.lpfilter     = lp_yn;
cfgpp.lpfreq       = lp_freq;
cfg.demean         = demean_yn;
cfg.baselinewindow = bsln_lim;
roi = ft_preprocessing(cfgpp, roi);

% Average ERPs
roi_erp = {};
n_trials = zeros([1 numel(cond_lab)]);
cfgavg = [];
cfgavg.keeptrials = 'yes';
for cond_ix = 1:numel(cond_lab)
    cfgavg.trials = find(cond_idx==cond_ix);
    roi_erp{cond_ix} = ft_timelockanalysis(cfgavg,roi);
    % Grab n_trials for design matrix
    n_trials(cond_ix) = size(roi_erp{cond_ix}.trial,1);
end

%% Contrast conditions
if strcmp(conditions,'DifOut') || numel(cond_lab)>2
    error(['Analysis not implemented for ' conditions ' yet, too many conditions: ' strjoin(cond_lab,',')]);
%     % Compute Win - Loss
%     cfgdif = [];
%     cfgdif.operation = 'subtract';
%     cfgdif.parameter = 'avg';
%     wn_ix = find(~cellfun(@isempty,strfind(cond_lab,'Wn')));
%     ls_ix = find(~cellfun(@isempty,strfind(cond_lab,'Ls')));
%     difwave = cell(size(wn_ix));
%     for dif_ix = 1:numel(wn_ix)
%         difwave{dif_ix} = roi_erp{wn_ix(dif_ix)};
%         % Compute difference in means (ERP difference wave)
%         difwave{dif_ix}.avg = difwave{dif_ix}.avg - roi_erp{ls_ix(dif_ix)}.avg;
%         difwave{dif_ix} = rmfield(difwave{dif_ix},'trial'); % remove trial since that doesn't make sense anymore
%         % Compute variance of difference in 2 RVs = var(X) + var(Y) - 2*covariance(X,Y)
%         wl_cov = cov(difwave{dif_ix}.avg, roi_erp{ls_ix(dif_ix)}.avg);
%         difwave{dif_ix}.var = difwave{dif_ix}.var + roi_erp{ls_ix(dif_ix)}.var - 2*wl_cov(1,2);
%     end
%     stat_lab = {strrep(cond_lab{wn_ix(1)},'Wn',''), strrep(cond_lab{wn_ix(2)},'Wn','')};
%     
%     cfgdw = [];
%     cfgdw.method = 'stats';
%     cfgdw.statistic = 'ttest2';
%     cfgdw.alpha = 0.05;
%     cfgdw.tail = 0;
%     cfgdw.parameter = 'avg';
%     cfgdw.design = [1 2];
%     stat = ft_timelockstatistics(cfgdw, difwave{:});
end

%% Run Statistics
% Create design matrix
% design = zeros(2,size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1));
% % Conditions (Independent Variable)
% design(1,1:size(roi_erp_con.trial,1)) = 1;
% design(1,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
% % Trial Numbers
% design(2,1:size(roi_erp_con.trial,1)) = 1;
% design(2,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
design = zeros(2,sum(n_trials));
for cond_ix = 1:numel(cond_lab)
    if cond_ix==1
        design(1,1:n_trials(cond_ix)) = cond_ix;                                % Conditions (Independent Variable)
        design(2,1:n_trials(cond_ix)) = 1:n_trials(cond_ix);                    % Trial Numbers
    else
        design(1,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= cond_ix; % Conditions (Independent Variable)
        design(2,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= 1:n_trials(cond_ix);
    end
end

% Prepare neighbors layout
% cfgn = [];
% cfgn.method  = 'distance';
% cfgn.layout  = 'ordered';
% cfgn.channel = elecs;
% neighbors    = ft_prepare_neighbours(cfgn,roi_erp_allch{1});
% for ch_ix = 1:numel(roi_erp{1}.label)
%     neighbors(ch_ix).label = roi_erp{1}.label{ch_ix};
%     neighbors(ch_ix).neighblabel = {};
% end

% Calculate statistics
cfg_stat.design           = design;
[stat] = ft_timelockstatistics(cfg_stat, roi_erp{:});

%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_filename);
save(data_out_filename,'-v7.3','roi_erp','stat');

% %% Plot Results
% stat_full = stat;
% roi_erp_full = roi_erp;
% 
% % Create a figure for each ROI
% for roi_ix = 1:numel(SBJ_vars.ch_lab.ROI)
%     % Select data to plot this ROI
%     cfgs = [];
%     cfgs.channel = SBJ_vars.ch_lab.ROI{roi_ix};
%     stat = ft_selectdata(cfgs,stat_full);
%     for an_ix = 1:numel(cond_lab)
%         roi_erp{an_ix} = ft_selectdata(cfgs,roi_erp_full{an_ix});
%     end
%     
%     % Plot parameters
%     roi_name = SBJ_vars.ch_lab.ROI{roi_ix}; if strcmp(roi_name(end),'*');roi_name=roi_name(1:end-1);end
%     fig_name = [SBJ '_ERP_stat_' conditions '_' roi_name '_' event_type];
%     [plot_rc,~] = fn_num_subplots(numel(stat.label));
%     if plot_rc(1)>1; fig_height=1; else fig_height=0.33; end;
%     
%     figure('Name',fig_name,'units','normalized',...
%         'outerposition',[0 0 1 fig_height],'Visible',fig_vis);
%     plot_info.fig        = gcf;
%     plot_info.x_step     = 0.25*roi.fsample;
%     plot_info.x_lab      = trial_lim_s(1):0.25:trial_lim_s(2);
%     plot_info.legend_loc = 'southeast';
%     plot_info.sig_alpha  = 0.2;
%     plot_info.sig_color  = [0.5 0.5 0.5];
%     % Stimulus plotting params
%     event_info.time      = -trial_lim_s(1)*roi.fsample;
%     event_info.name      = 'stim';
%     event_info.width     = 2;
%     event_info.color     = 'k';
%     event_info.style     = '--';
%     % Condition plotting params
%     cond_info.name       = cond_lab;
%     cond_info.style      = cond_style;
%     cond_info.color      = cond_colors;
%     cond_info.alpha      = [0.5 0.5];
%     
%     % Plot each channel within this ROI
%     for ch_ix = 1:numel(stat.label)
%         subplot(plot_rc(1),plot_rc(2),ch_ix);
%         plot_info.ax         = gca;
%         plot_info.title      = stat.label{ch_ix};
%         if ch_ix==1; plot_info.legend=1; else plot_info.legend=0; end;
%         
%         % Compute means and variance
%         means = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
%         var = NaN([numel(cond_lab) size(roi_erp{1}.avg,2)]);
%         for an_ix = 1:numel(cond_lab)
%             means(an_ix,:) = roi_erp{an_ix}.avg(ch_ix,:);
%             var(an_ix,:) = squeeze(std(roi_erp{an_ix}.trial(:,ch_ix,:),[],1)./sqrt(size(roi_erp{an_ix}.trial,1)))';
%         end
%         % Find significant time periods
%         if sum(stat.mask(ch_ix,:))>0
%             mask_chunks = fn_find_chunks(stat.mask(ch_ix,:));
%             sig_chunks = mask_chunks;
%             sig_chunks(stat.mask(ch_ix,sig_chunks(:,1))==0,:) = [];
%             % If stat and roi_erp aren't on same time axis, adjust sig_chunk indices
%             if (size(stat.time,2)~=size(roi_erp{1}.time,2)) || (sum(stat.time==roi_erp{1}.time)~=numel(stat.time))
%                 for chunk_ix = 1:size(sig_chunks,1)
%                     sig_chunks(chunk_ix,1) = find(roi_erp{1}.time==stat.time(sig_chunks(chunk_ix,1)));
%                     sig_chunks(chunk_ix,2) = find(roi_erp{1}.time==stat.time(sig_chunks(chunk_ix,2)));
%                 end
%             end
%             fprintf('%i SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',size(sig_chunks,1));
%             fn_plot_ts_error_bar_sig(plot_info,means,var,sig_chunks,event_info,cond_info);
%         else
%             fprintf('NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n');
%             fn_plot_ts_error_bar(plot_info,means,var,event_info,cond_info);
%         end
%     end
%     clear stat roi_erp
%     
%     % Save figure
%     if save_fig
%         fig_dir = ['/home/knight/hoycw/PRJ_Error/results/ERP/' SBJ '/' conditions '/'];
%         if ~exist(fig_dir,'dir')
%             mkdir(fig_dir);
%         end
%         fig_filename = [fig_dir fig_name '.' fig_filetype];
%         fprintf('Saving %s\n',fig_filename);
%         saveas(gcf,fig_filename);
%         %eval(['export_fig ' fig_filename]);
%     end
% end
%%
% % Plot ERPs
% roi_erp_con.mask = stat.mask;
% roi_erp_inc.mask = stat.mask;
% cfgp = [];
% cfgp.showlabels = 'yes';
% cfgp.parameter = 'avg';
% % cfgp.layout = 'ordered';
% cfgp.maskparameter = 'mask';
% ft_singleplotER(cfgp, roi_erp_con, roi_erp_inc);%stat); %roi_erp_con, roi_erp_inc, 

end
