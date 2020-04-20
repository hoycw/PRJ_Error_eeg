function SBJ04c_ERP_grp_stats_LME_P2P(SBJ_id,proc_id,an_id,stat_id,varargin)
% Compute linear mixed-effects model on SBJ average data
%   Should only be used for peak-to-peak FRN
%   Jackknife DF correction unclear, so single SBJ with no peak = NaN
%   Only runs for one channel right now...
% INPUTS:
% OUTPUTS:

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
        elseif strcmp(varargin{v},'plot_erps')
            plot_erps = varargin{v+1};
        elseif strcmp(varargin{v},'plot_peaks')
            plot_peaks = varargin{v+1};
        else
            error(['Unknown varargin ' num2str(v) ': ' varargin{v}]);
        end
    end
end

% Define default options
if ~exist('fig_vis','var');    fig_vis = 'on'; end
if ~exist('fig_ftype','var');  fig_ftype = 'png'; end
if ~exist('save_fig','var');   save_fig = 1; end
if ~exist('plot_erps','var');  plot_erps = 0; end
if ~exist('plot_peaks','var'); plot_peaks = 1; end
fig_dir = [root_dir 'PRJ_Error_eeg/results/ERP/' an_id '/' stat_id '/GLM_inputs/'];
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
if ~strcmp(st.measure,'p2p') || ~strcmp(st.an_style,'lme'); error('This LME should only be run on p2p!'); end
if plot_peaks && ~strcmp(st.measure,'p2p'); error('why plot peaks ifnot using p2p?'); end
if strcmp(st.grp_method,'jackknife'); error('jackknife not figrued out for GLM yet!'); end

% Select SBJs
SBJs = load_SBJ_file(SBJ_id);

model_id = [st.model_lab '_' st.trial_cond{1}];
[reg_lab, ~, ~, ~]     = fn_regressor_label_styles(st.model_lab);
[cond_lab, ~, cond_colors, ~, ~] = fn_condition_label_styles(st.trial_cond{1});

%% Load Behavior
bhvs          = cell(size(SBJs));
full_cond_idx = cell(size(SBJs));
n_trials      = zeros([numel(SBJs) 1]);
for s = 1:numel(SBJs)
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
cfgs  = []; cfgs.latency = st.stat_lim;
model = zeros([numel(cond_lab)*numel(SBJs) numel(reg_lab)]);
sbj_factor = [];
full_model_ix = 0;
for s = 1:numel(SBJs)
    % Load data
    fprintf('========================== Processing %s ==========================\n',SBJs{s});
    load([root_dir 'PRJ_Error_eeg/data/',SBJs{s},'/04_proc/',SBJs{s},'_',an_id,'.mat'],'roi');
    if numel(roi.label)>1; error('only ready for one channel in P2P GLM!'); end
    
    % Select time and trials of interest
    cfgs.trials  = find(full_cond_idx{s});
    st_roi = ft_selectdata(cfgs, roi);
    
    if s==1
        % Initialize matrices now that we know time axis
        time_vec = st_roi.time{1};
        ch_list  = st_roi.label;
        erps  = nan([numel(cond_lab) numel(SBJs) numel(ch_list) numel(time_vec)]);
    end
    
    % Load RL Model
    tmp = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/04_proc/' SBJs{s} '_model_' model_id '.mat']);
    
    % Z-score SBJ model regressors
    sbj_model = NaN(size(tmp.model));
    if st.z_reg
        for reg_ix = 1:numel(reg_lab)
            sbj_model(:,reg_ix) = ...
                (tmp.model(:,reg_ix)-nanmean(tmp.model(:,reg_ix)))./nanstd(tmp.model(:,reg_ix));
        end
    else
        sbj_model = tmp.model;
    end
    
    % Compute ERPs and average model within condition
    cond_idx = fn_condition_index(cond_lab, bhvs{s});
    for cond_ix = 1:numel(cond_lab)
        cond_trial_ix = find(cond_idx==cond_ix);
        % Compute ERP
        trials = nan([numel(ch_list) numel(cond_trial_ix) numel(time_vec)]);
        for trl_ix = 1:numel(cond_trial_ix)
            trials(:,trl_ix,:) = st_roi.trial{cond_trial_ix(trl_ix)};
        end
        erps(cond_ix,s,:,:) = mean(trials,2);
        
        % Average model
        full_model_ix = full_model_ix + 1;
        model(full_model_ix, :) = mean(sbj_model(cond_trial_ix,:),1);
    end
    
    % Track SBJ
    sbj_factor = [sbj_factor; s*ones([numel(cond_lab) 1])];
    
    clear tmp roi st_roi sbj_model trials cond_idx cond_trial_ix
end

%% Compute Peak-to-Peak Data
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
    % Find Peak Search Range
    pk_rng = zeros([2 2]);
    pk_sign_str = cell([2 1]);
    for pk_ix = 1:2
        for lim_ix = 1:2
            [~, pk_rng(pk_ix,lim_ix)] = min(abs(time_vec-st.pk_lim(pk_ix,lim_ix)));
        end
        if st.pk_sign(pk_ix)==1
            pk_sign_str{pk_ix} = 'positive';
        elseif st.pk_sign(pk_ix)==-1
            pk_sign_str{pk_ix} = 'negative';
        else
            error('st.pk_sign not 1/-1!');
        end
    end
    % Detect peaks and compute amplitude difference
    if strcmp(st.grp_method,'jackknife')
        data = zeros([numel(cond_lab) numel(SBJs) numel(ch_list)]);
    elseif strcmp(st.grp_method,'sbj')
        data = nan([numel(cond_lab) numel(SBJs) numel(ch_list)]);
    end
    pk_amp   = nan([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
    pk_times = nan([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
    mult_pks = false([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
    miss_pks = false([numel(cond_lab) numel(SBJs) numel(ch_list) 2]);
    bad_pks  = false([numel(cond_lab) numel(SBJs) numel(ch_list)]);
    for s = 1:numel(SBJs)
        for cond_ix = 1:numel(cond_lab)
            for ch_ix = 1:numel(ch_list)
                % Find positive and negative peak amplitudes and latencies
                for pk_ix = 1:2
                    [tmp_amp, tmp_lat] = findpeaks(st.pk_sign(pk_ix)*...
                        squeeze(gavg(cond_ix,s,ch_ix,pk_rng(pk_ix,1):pk_rng(pk_ix,2))));
                    
                    % Select correct peak
                    if ~isempty(tmp_amp)
                        % Remove possible 2nd peaks preceding 1st peak
                        if pk_ix==2 && any(time_vec(tmp_lat+pk_rng(pk_ix,1)-1)<pk_times(cond_ix,s,ch_ix,1))
                            tmp_amp(time_vec(tmp_lat+pk_rng(pk_ix,1)-1)<pk_times(cond_ix,s,ch_ix,1)) = [];
                            tmp_lat(time_vec(tmp_lat+pk_rng(pk_ix,1)-1)<pk_times(cond_ix,s,ch_ix,1)) = [];
                        end
                        if numel(tmp_amp)>1
                            % Multiple peaks detected, take maximum amplitude
                            %   + peaks: most positive
                            %   - peaks: most positive of negatives, so lowest non-inverted amplitude
                            [~, amp_ix] = max(tmp_amp);
                            mult_pks(cond_ix,s,ch_ix,pk_ix) = true;
                            fprintf('\tMultiple %s peaks detected for %s %s!\n',...
                                pk_sign_str{pk_ix},SBJs{s},cond_lab{cond_ix});
                        else
                            amp_ix = 1;
                        end
                        pk_amp(cond_ix,s,ch_ix,pk_ix) = st.pk_sign(pk_ix)*tmp_amp(amp_ix);
                        pk_times(cond_ix,s,ch_ix,pk_ix) = time_vec(tmp_lat(amp_ix)+pk_rng(pk_ix,1)-1);
                    else
                        % No peak detected
                        miss_pks(cond_ix,s,ch_ix,pk_ix) = true;
                        fprintf(2,'\tNo %s peak detected for %s %s!\n',...
                            pk_sign_str{pk_ix},SBJs{s},cond_lab{cond_ix});
                    end
                end
                
                % Compute peak-to-peak amplitude difference
                if ~miss_pks(cond_ix,s,ch_ix)
                    amp_diff = pk_amp(cond_ix,s,ch_ix,1)-pk_amp(cond_ix,s,ch_ix,2);
                    % Toss differences that go in the wrong direction or order
                    %   This indicates a complex, multi-peak waveform...
                    if amp_diff > 0 && pk_times(cond_ix,s,ch_ix,1)<pk_times(cond_ix,s,ch_ix,2)
                        data(cond_ix,s,ch_ix) = amp_diff;
                    else
                        bad_pks(cond_ix,s,ch_ix) = true;
                        fprintf(2,'\tBad %s peak detected for %s %s!\n',...
                            pk_sign_str{pk_ix},SBJs{s},cond_lab{cond_ix});
                    end
                end
                
                % Plot Peak Detection Errors
                if plot_peaks && (any(mult_pks(cond_ix,s,ch_ix,:)) || ...
                        any(miss_pks(cond_ix,s,ch_ix,:)) || bad_pks(cond_ix,s,ch_ix))
                    fig_name = ['bad_peaks_' st.grp_method '_' SBJs{s} ...
                        '_' cond_lab{cond_ix} '_' ch_list{ch_ix}];
                    figure('Name',fig_name,'Visible',fig_vis);
                    plot(time_vec(pk_rng(1,1):pk_rng(2,2)),...
                        squeeze(gavg(cond_ix,s,ch_ix,pk_rng(1,1):pk_rng(2,2))),...
                        'Color',cond_colors{cond_ix});
                    hold on;
                    for pk_ix = 1:2
                        [tmp_amp, tmp_lat] = ...
                            findpeaks(st.pk_sign(pk_ix)*...
                            squeeze(gavg(cond_ix,s,ch_ix,pk_rng(pk_ix,1):pk_rng(pk_ix,2))));
                        if st.pk_sign(pk_ix)==1; mrkr_color = 'b'; else mrkr_color = 'r'; end
                        if ~miss_pks(cond_ix,s,ch_ix,pk_ix)
                            scatter(time_vec(tmp_lat+pk_rng(pk_ix,1)-1),...
                                st.pk_sign(pk_ix)*tmp_amp, 'o', mrkr_color);
                        end
                    end
                    bad_str = [];
                    if any(mult_pks(cond_ix,s,ch_ix,:))
                        bad_str = [bad_str ' mult' num2str(sum(mult_pks(cond_ix,s,ch_ix,:))) ';'];
                    end
                    if bad_pks(cond_ix,s,ch_ix)
                        bad_str = [bad_str ' bad;'];
                    end
                    if any(miss_pks(cond_ix,s,ch_ix,:))
                        bad_str = [bad_str ' miss' num2str(sum(miss_pks(cond_ix,s,ch_ix,:))) ';'];
                    else
                        line([pk_times(cond_ix,s,ch_ix,1) pk_times(cond_ix,s,ch_ix,2)],...
                            [pk_amp(cond_ix,s,ch_ix,1) pk_amp(cond_ix,s,ch_ix,2)],...
                            'Color','k');
                    end
                    title([SBJs{s} ' ' cond_lab{cond_ix} ': ' bad_str]);
                    xlabel('Time (s)'); ylabel('Amplitude (uV)');
                    set(gca,'FontSize',16);
                    hold off;
                    if save_fig
                        bad_fig_dir = [fig_dir 'bad_peak_detection/'];
                        if ~exist(bad_fig_dir,'dir'); mkdir(bad_fig_dir); end
                        fig_fname = [bad_fig_dir fig_name '.' fig_ftype];
                        %fprintf('Saving %s\n',fig_fname);
                        saveas(gcf,fig_fname);
                    else
                        pause;
                    end
                end
            end
        end
    end
else
    error(['unknown st.method for ' stat_id]);
end

% Reshape means for GLM
st_data = reshape(data,[numel(cond_lab)*numel(SBJs), numel(ch_list)]);
if any(isnan(st_data(:)))
    if strcmp(st.grp_method,'jackknife')
        error([num2str(sum(isnan(st_data(:)))) ' NaNs in GLM data! Cannot correct F stat in imbalanced cells!']);
    else
        warning([num2str(sum(isnan(st_data(:)))) ' NaNs in GLM data!']);
    end
end

%% Plot Data
% Plot ERPs (original or jackknife)
if plot_erps
    for ch_ix = 1:numel(ch_list)
        fig_name = ['ERPs_' st.grp_method '_' ch_list{ch_ix} '_' SBJ_id];
        figure('Name',fig_name,'Visible',fig_vis);
        tmp = gavg(:,:,ch_ix,:);
        ylims = [min(tmp(:)) max(tmp(:))];
        for cond_ix = 1:numel(cond_lab)
%             subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix); hold on;
            subplot(2,3,cond_ix); hold on;
            title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
            plot(time_vec,squeeze(gavg(cond_ix,:,ch_ix,:)),...
                'Color',cond_colors{cond_ix});%,'LineStyle',cond_styles{cond_ix});
            ylim(ylims);
            fn_min_white_space(gca);
        end
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
    end
end

% Plot peaks in peak range
if plot_peaks
    for ch_ix = 1:numel(ch_list)
        fig_name = ['ERPs_peaks_' st.grp_method '_'  ch_list{ch_ix} '_' SBJ_id];
        figure('Name',fig_name,'Visible',fig_vis);
        tmp = pk_amp(:,:,ch_ix,:);
        ylims = [min(tmp(:)) max(tmp(:))];
        tmp = pk_times(:,:,ch_ix,:);
        xlims = [min(tmp(:)) max(tmp(:))];
        for cond_ix = 1:numel(cond_lab)
%             subplot(numel(grp_cond_lab{1}),numel(grp_cond_lab{2}),cond_ix); hold on;
            subplot(2,3,cond_ix); hold on;
            title([ch_list{ch_ix} ': ' cond_lab{cond_ix}]);
            plot(time_vec(pk_rng(1,1):pk_rng(2,2)),...
                squeeze(gavg(cond_ix,:,ch_ix,pk_rng(1,1):pk_rng(2,2))),...
                'Color',cond_colors{cond_ix});%,'LineStyle',cond_styles{cond_ix});
            scatter(pk_times(cond_ix,:,ch_ix,1), pk_amp(cond_ix,:,ch_ix,1), 'o', 'b');
            scatter(pk_times(cond_ix,:,ch_ix,2), pk_amp(cond_ix,:,ch_ix,2), 'o', 'r');
            for s = 1:numel(SBJs)
                line([pk_times(cond_ix,s,ch_ix,1) pk_times(cond_ix,s,ch_ix,2)],...
                    [pk_amp(cond_ix,s,ch_ix,1) pk_amp(cond_ix,s,ch_ix,2)],...
                    'Color','k');
            end
            ylim(ylims);
            xlim(xlims);
            fn_min_white_space(gca);
        end
        if save_fig
            fig_fname = [fig_dir fig_name '.' fig_ftype];
            fprintf('Saving %s\n',fig_fname);
            saveas(gcf,fig_fname);
        end
    end
end

%% Compute Statistics
fprintf('========================== Running Stats ==========================\n');
tic
% Build Model Table
tbl = table;
for reg_ix = 1:numel(reg_lab)
    tbl.(reg_lab{reg_ix}) = model(:,reg_ix);
end
tbl.SBJ = categorical(sbj_factor);

% Create Model Formula
reg_formula = strjoin(reg_lab,' + ');
formula = ['ERP ~ ' reg_formula ' + (1|SBJ)'];

% Run Model
if strcmp(st.measure,'ts')
%     glm = cell(size(time_vec));
%     pvals = nan([numel(reg_lab) numel(time_vec)]);
%     for t_ix = 1:numel(time_vec)
%         tbl.ERP = data(:,t_ix);
%         glm{t_ix} = fitglm(tbl,formula);
%         pvals(:,t_ix) = glm{t_ix}.Coefficients.pValue(2:end);
%     end
elseif any(strcmp(st.measure,{'mean','p2p'}))
    lme   = cell(size(ch_list));
    pvals = nan([numel(reg_lab) numel(ch_list)]);
    for ch_ix = 1:numel(ch_list)
        tbl.ERP    = st_data(:,ch_ix);
        lme{ch_ix} = fitlme(tbl,formula);
        pvals(:,ch_ix) = lme{ch_ix}.Coefficients.pValue(2:end);
        
        % Adjust p for jackknife
        if strcmp(st.grp_method,'jackknife')
            % This adjustment doesn't work! fitglm does not return an F
            % stat, and it's unclear what DF should be used to correct anyways
%             for reg_ix = 1:numel(reg_lab)
%                 % Adjusted F = F / (n-1)^2;
%                 pvals(reg_ix,ch_ix) = 1-fcdf(glm{1}.anova.FStat(reg_ix+1)/(numel(SBJs)-1)^2,...
%                     glm{1}.anova.DF1(reg_ix+1), glm{1}.anova.DF2(reg_ix+1));
%             end
        end
    end
else
    error(['Unknown st.measure: ' st.measure]);
end

% Correct for Multiple Comparisons
if strcmp(st.mcp_method,'FDR')
    [~, ~, ~, qvals] = fdr_bh(reshape(pvals,[size(pvals,1)*size(pvals,2) 1]));
    qvals = reshape(qvals,[size(pvals,1) size(pvals,2)]);
else
    error(['Unknown method for multiple comparison correction: ' st.mcp_method]);
end

fprintf('\t\t Stats Complete:');
toc

%% Save Results
stat_out_dir = [root_dir 'PRJ_Error_eeg/data/GRP/'];
if ~exist(stat_out_dir,'dir')
    mkdir(stat_out_dir);
end
stat_out_fname = [stat_out_dir SBJ_id '_' stat_id '_' an_id '.mat'];
fprintf('Saving %s\n',stat_out_fname);
save(stat_out_fname,'-v7.3','lme','qvals','SBJs','time_vec','ch_list','pk_amp','pk_times','data');

end
