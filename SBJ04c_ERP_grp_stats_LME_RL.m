function SBJ04c_ERP_grp_stats_LME_RL(SBJs,proc_id,an_id,stat_id)
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

[cond_lab, ~, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Behavior, Compute Model
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
pWin          = cell(size(SBJs));
sPE           = cell(size(SBJs));
uPE           = cell(size(SBJs));
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
    
    % Run Logistic Regression for Win Prediction
    % Select Data (fit on everything except surprise since no outcome)
    s_idx = fn_condition_index({'Su'},bhvs{s});
    X = bhvs{s}.tol(~s_idx);
    y = double(bhvs{s}.hit(~s_idx));
    
    % Logistic regression
    betas = glmfit(X,y,'binomial','link','logit');
    
    z = betas(1) + (bhvs{s}.tol * betas(2));
    pWin{s} = 1 ./ (1+exp(-z));
    expected_score = pWin{s}*2 - 1;
    sPE{s} = double(bhvs{s}.score)/100 - expected_score;
    uPE{s} = abs(sPE{s});
    clear tmp
end

%% Load Data and Build Model
full_trl_ix = 0;
cfgs = []; cfgs.latency = st.stat_lim;
%model = nan([sum(n_trials) numel(st.factors)+1]);
sbj_factor = zeros([sum(n_trials) 1]);
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
        if strcmp(st.measure,'ts')
            data  = nan([sum(n_trials) numel(time_vec)]);
        elseif strcmp(st.measure,'mean')
            data  = nan(size(n_trials));
        else; error(['unknown st.measure: ' st.measure]);
        end
        
        % Track SBJ
        sbj_factor(1:n_trials(s),end) = s*ones([n_trials(s) 1]);
    else
        % Track SBJ
        sbj_factor(sum(n_trials(1:s-1))+1:sum(n_trials(1:s))) = s*ones([n_trials(s) 1]);
    end
    
    % Add data and model
    for trl_ix = 1:numel(st_roi.trial)
        full_trl_ix = full_trl_ix + 1;
        if strcmp(st.measure,'ts')
            data(full_trl_ix,:) = st_roi.trial{trl_ix};
        elseif strcmp(st.measure,'mean')
            data(full_trl_ix) = mean(st_roi.trial{trl_ix});
        else; error(['unknown st.measure: ' st.measure]);
        end
    end
    
    clear roi st_roi
end

%% Build table
fprintf('========================== Running Stats ==========================\n');
tic
if strcmp(st.measure,'ts')
    lme = cell(size(time_vec));
    pvals = nan([numel(st.factors) numel(time_vec)]);
    for t_ix = 1:numel(time_vec)
        tbl = table(data(:,t_ix),vertcat(pWin{:}),vertcat(sPE{:}),vertcat(uPE{:}),sbj_factor,...
            'VariableNames',{'ERP',st.factors{:},'SBJ'});
        %     for grp_ix = 1:numel(st.factors)
        %         if st.categorical(grp_ix)
        %             tbl.(st.factors{grp_ix}) = categorical(tbl.(st.factors{grp_ix}));
        %         end
        %     end
        tbl.SBJ = categorical(tbl.SBJ);
        lme{t_ix} = fitlme(tbl,st.formula);
        for fct_ix = 1:numel(st.factors)
            pvals(fct_ix,t_ix) = lme{t_ix}.Coefficients.pValue(fct_ix+1);
        end
    end
elseif strcmp(st.measure,'mean')
    lme = {};
    pvals = nan([numel(st.factors) 1]);
    tbl = table(data,vertcat(pWin{:}),vertcat(sPE{:}),vertcat(uPE{:}),sbj_factor,...
        'VariableNames',{'ERP',st.factors{:},'SBJ'});
    tbl.SBJ = categorical(tbl.SBJ);
    lme{1} = fitlme(tbl,st.formula);
    for fct_ix = 1:numel(st.factors)
        pvals(fct_ix) = lme{1}.Coefficients.pValue(fct_ix+1);
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
