function SBJ04c_ERP_grp_stats_LME(SBJs,proc_id,an_id,stat_id)
% Run Mixed-Effects Linear model on all SBJ and trials
% INPUTS:
%   SBJs [cell array] - ID list of subjects to run
%   proc_id [str] - ID of preprocessing pipeline
%   an_id [str] - ID of the analysis parameters to use
%   stat_id [str] - ID of the stats parameters to use
% OUTPUTS:
%   lme [cell array] - LinearMixedModel output class for each time point

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
elseif exist('Users/aasthashah/', 'dir'); root_dir = 'Users/aasthashah/Desktop/'; ft_dir = 'Users/aasthashah/Applications/fieldtrip';
else; root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);
addpath([app_dir 'fieldtrip/']);
ft_defaults

%% Load Data 
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Select Conditions of Interest
[grp_lab, ~, ~]  = fn_group_label_styles(st.model_lab);
[cond_lab, ~, ~, ~] = fn_condition_label_styles(st.model_lab);
fact_cond_lab = cell(size(st.factors));
for fct_ix = 1:numel(st.factors)
    [fact_cond_lab{fct_ix}, ~, ~, ~] = fn_condition_label_styles(st.factors{fct_ix});
    if numel(fact_cond_lab{fct_ix}) > 2; error('coding not ready for 3 conditions'); end
end

%% Load Behavior
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros(size(SBJs));
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{1} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load data
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' ...
        SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    bhvs{s} = tmp.bhv;
    
    % Select Conditions of Interest
    full_cond_idx{s} = fn_condition_index(cond_lab, bhvs{s});
    bhv_fields = fieldnames(bhvs{s});
    orig_n_trials = numel(bhvs{s}.trl_n);
    for f_ix = 1:numel(bhv_fields)
        if numel(bhvs{s}.(bhv_fields{f_ix}))==orig_n_trials
            bhvs{s}.(bhv_fields{f_ix}) = bhvs{s}.(bhv_fields{f_ix})(full_cond_idx{s}~=0);
        end
    end
    
    n_trials(s) = numel(bhvs{s}.trl_n);
    clear tmp
end

%% Load Data and Build Model
full_trl_ix = 0;
cfgs = []; cfgs.latency = st.stat_lim;
model = nan([sum(n_trials) numel(st.factors)+1]);
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',an_id,'.mat'],'roi');
    if numel(roi.label)>1; error('assuming single channel for now!'); end
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_roi = ft_selectdata(cfgs, roi);
    
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = st_roi.time{1};
        data  = nan([sum(n_trials) numel(time_vec)]);
        
        % Create model
        for fct_ix = 1:numel(fact_cond_lab)
            % 1 for first condition, leave 0 for second
            model(1:n_trials(s),fct_ix) = fn_condition_index(fact_cond_lab{fct_ix}, bhvs{s});
        end
        model(1:n_trials(s),end) = s*ones([n_trials(s) 1]);
    else
        % Create model
        for fct_ix = 1:numel(fact_cond_lab)
            model(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),fct_ix) = fn_condition_index(fact_cond_lab{fct_ix}, bhvs{s});
        end
        model(sum(n_trials(1:s-1))+1:sum(n_trials(1:s)),end) = s*ones([n_trials(s) 1]);
    end
    
    % Add data and model
    for trl_ix = 1:numel(st_roi.trial)
        full_trl_ix = full_trl_ix + 1;
        data(full_trl_ix,:) = st_roi.trial{trl_ix};
    end
    
    clear roi st_roi
end

%% Build table
fprintf('========================== Running Stats ==========================\n');
tic
lme = cell(size(time_vec));
pvals = nan([numel(grp_lab) numel(time_vec)]);
for t_ix = 1:numel(time_vec)
    tbl = table(data(:,t_ix),model(:,1),model(:,2),model(:,3),'VariableNames',{'ERP',st.factors{:},'SBJ'});
    lme{t_ix} = fitlme(tbl,st.formula);
    for grp_ix = 1:numel(grp_lab)
        label_ix = strcmp(lme{t_ix}.CoefficientNames,strrep(grp_lab{grp_ix},'*',':'));
        pvals(grp_ix,t_ix) = lme{t_ix}.Coefficients.pValue(label_ix);
    end
end
fprintf('\t\t Stats Complete:');
toc

%% Save Results
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(stat_out_dir,'dir')
    [~] = mkdir(stat_out_dir);
end
stat_out_fname = [stat_out_dir 'GRP_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','SBJs');

end
