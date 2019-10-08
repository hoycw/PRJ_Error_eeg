function SBJ03c_ERP_grp_stats(SBJs,conditions,proc_id,an_id)
% Compute grand average group ERP from SBJ ERPs:
%   Re-align data to event, select channels and epoch, filter, average, run stats, save
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
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/', ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data 
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);

% Select Conditions of Interest
[cond_lab, ~, ~, ~] = fn_condition_label_styles(conditions);

% Load Data
SBJs_vars = cell(size(SBJs));
rois      = cell(size(SBJs));
w2s       = cell(size(SBJs));
for s = 1:length(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    SBJs_vars{s} = SBJ_vars;
    
    tmp = load([SBJs_vars{s}.dirs.SBJ,'04_proc/',SBJs{s},'_',an_id,'.mat']);
    rois{s} = tmp.roi; w2s{s} = tmp.w2;
    clear SBJ_vars
end

%% Compute Grand Average Group ERPs
st_data = nan([numel(SBJs) numel(w2{1}.label) numel(w2{1}.time)]);
for s = 1:numel(SBJs)
    st_data(s,:,:) = w2{s}.zscore;
end
if any(isnan(st_data(:))); error('NaN in ANOVA data!'); end

% grp_erp = cell(size(cond_lab));
% for cond_ix = 1:numel(cond_lab)
%     grp_erp{cond_ix} = ft_timelockgrandaverage(cfg_gavg,roi_erps{:,cond_ix});
% end

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
% % Create design matrix
% design = zeros(2,numel(SBJs));
% for cond_ix = 1:numel(cond_lab)
%     design(1,cond_ix
%     if cond_ix==1
%         design(1,1:numel(SBJs)) = cond_ix;          % Conditions (Independent Variable)
%         design(2,1:numel(SBJs)) = 1:numel(SBJs);	% SBJ Numbers
%     else
%         design(1,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= cond_ix; % Conditions (Independent Variable)
%         design(2,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= 1:n_trials(cond_ix);
%     end
% end
% 
% % Prepare neighbors layout
% % cfgn = [];
% % cfgn.method  = 'distance';
% % cfgn.layout  = 'ordered';
% % cfgn.channel = elecs;
% % neighbors    = ft_prepare_neighbours(cfgn,roi_erp_allch{1});
% % for ch_ix = 1:numel(roi_erp{1}.label)
% %     neighbors(ch_ix).label = roi_erp{1}.label{ch_ix};
% %     neighbors(ch_ix).neighblabel = {};
% % end
% 
% % Calculate statistics
% cfg_stat.design = design;
% % [stat] = ft_timelockstatistics(cfg_stat, roi_erp{:});
% warning('HACK!!!! MEX file work around, just copying roi_erp instead of running real stats!!!');
% stat = roi_erp{1};
% stat.mask = zeros([numel(stat.label) numel(stat.time)]);

%% Save Results
data_out_fname = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_',an_id,'.mat');
fprintf('Saving %s\n',data_out_fname);
save(data_out_fname,'-v7.3','grp_erp','SBJs');

end
