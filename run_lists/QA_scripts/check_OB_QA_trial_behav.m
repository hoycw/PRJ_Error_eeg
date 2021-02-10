%% Check trial counts for Oddball Task
% Task parameters
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

proc_id    = 'odd_full_ft';
SBJ_id     = 'goodEEG';
conditions = 'OB';

outlier_std_thresh = 3;

% Plotting colors
%   1 = good, 2 = suspect, 3 = likely toss, 4 = bad
colors = [...
    0 0 0; ...good- black
    1 0.5 0; ...suspect- orange
    1 0 1; ...likely toss - magenta
    1 0 0]; % bad- red
sz = 75;

%% Load data
% ICs tossed, bad trials tossed, channels tossed, bad epochs, bad RTs
SBJs = fn_load_SBJ_list(SBJ_id);
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
[cond_lab,~,~,~,~] = fn_condition_label_styles(conditions);

n_trials   = nan([numel(SBJs) numel(cond_lab)]);
n_miss_hit = zeros([numel(SBJs) numel(cond_lab)]);
n_false_hit= zeros([numel(SBJs) numel(cond_lab)]);
n_rt_out   = nan([numel(SBJs) numel(cond_lab)]);
n_bad_ics  = nan(size(SBJs));
n_bad_trl  = nan(size(SBJs));
accuracy   = nan(size(SBJs));
SBJ_colors = zeros([numel(SBJs) 3]);
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    n_bad_ics(s) = numel(SBJ_vars.ica_reject);
    n_bad_trl(s) = numel(SBJ_vars.trial_reject_ix_oddball);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    
%     fnames = fieldnames(bhv);
%     for f = 1:numel(fnames)
%         if ~strcmp(fnames{f},'blk')
%             bhv.(fnames{f}) = bhv.(fnames{f})(bhv.blk~=1);
%         end
%     end
%     bhv.blk = bhv.blk(bhv.blk~=1);
    
    % Behavioral Metrics
    cond_idx = fn_condition_index(cond_lab,bhv);
    for cond_ix = 1:numel(cond_lab)
        n_trials(s,cond_ix) = sum(cond_idx==cond_ix);
        if strcmp(cond_lab{cond_ix},'Tar')
            n_miss_hit(s,cond_ix)   = sum(bhv.hit(cond_idx==cond_ix)~=1);
        else
            n_false_hit(s,cond_ix)   = sum(bhv.miss(cond_idx==cond_ix));
        end
        rts = bhv.rt(cond_idx==cond_ix & bhv.rt~=-1);
        n_rt_out(s,cond_ix) = sum(rts<proc.rt_bounds(1)) + sum(rts>proc.rt_bounds(2));
    end
    accuracy(s) = 1 - sum(n_miss_hit(s,:)+n_false_hit(s,:))/sum(n_trials(s,:));
    
    clear SBJ_vars bhv cond_idx
end

% Compute outliers for bad trial count and accuracy
avg_bad_trl = nanmean(n_bad_trl);
std_bad_trl = nanstd(n_bad_trl);
avg_acc     = nanmean(accuracy);
std_acc     = nanstd(accuracy);
for s = 1:numel(SBJs)
    if accuracy(s) <= avg_acc - std_acc*outlier_std_thresh
        fprintf(2,'%s bad for accuracy %.3f\n',SBJs{s},accuracy(s));
        SBJ_colors(s,:) = [1 0 0];
    end
    if n_bad_trl(s) >= avg_bad_trl + std_bad_trl*outlier_std_thresh
        fprintf(2,'%s bad for n_bad_trl %d\n',SBJs{s},n_bad_trl(s));
        SBJ_colors(s,:) = [1 0 0];
    end
end

%% Data Qualty: Trial count, bad ICs, bad trials
figure;

% Bad Trials
subplot(4,numel(cond_lab),1); hold on;
scatter(1:numel(SBJs),n_bad_trl,sz,SBJ_colors,'filled');
%     line(xlim,[2 2],'Color','r');
ylims = ylim;
set(gca,'YLim',[0 ylims(2)]);
set(gca,'XLim',[0 numel(SBJs)+1]);
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
ylabel('# Bad Trials');
title('Bad Trial Count');
set(gca,'FontSize',14);

% Bad ICs
subplot(4,numel(cond_lab),2); hold on;
scatter(1:numel(SBJs),n_bad_ics,sz,SBJ_colors,'filled');
%     line(xlim,[30 30],'Color','r');
ylims = ylim;
set(gca,'YLim',[0 ylims(2)]);
set(gca,'XLim',[0 numel(SBJs)+1]);
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
ylabel('# Bad ICs');
title('Bad IC Count');
set(gca,'FontSize',14);

% Total Inaccurate Trials
subplot(4,numel(cond_lab),3); hold on;
scatter(1:numel(SBJs),accuracy,sz,SBJ_colors,'filled');
set(gca,'YLim',[0.95 1]);
set(gca,'XLim',[0 numel(SBJs)+1]);
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
ylabel('Accuracy');
title('Accuracy');
set(gca,'FontSize',14);

for cond_ix = 1:numel(cond_lab)
    % Trial Count
    subplot(4,numel(cond_lab),3+cond_ix); hold on;
    scatter(1:numel(SBJs),n_trials(:,cond_ix),sz,SBJ_colors,'filled');
    line(xlim,[30 30],'Color','r');
    ylims = ylim;
    set(gca,'YLim',[0 ylims(2)]);
    set(gca,'XLim',[0 numel(SBJs)+1]);
    set(gca,'XTick',1:numel(SBJs));
    set(gca,'XTickLabel',SBJs);
    xtickangle(45);
    xlabel('SBJ');
    ylabel('# Trials');
    title(cond_lab{cond_ix});
    set(gca,'FontSize',14);
    
    % Bad RTs
    subplot(4,numel(cond_lab),6+cond_ix); hold on;
    scatter(1:numel(SBJs),n_rt_out(:,cond_ix),sz,SBJ_colors,'filled');
    line(xlim,[2 2],'Color','r');
    set(gca,'XLim',[0 numel(SBJs)+1]);
    set(gca,'XTick',1:numel(SBJs));
    set(gca,'XTickLabel',SBJs);
    xtickangle(45);
    xlabel('SBJ');
    ylabel('# Bad RTs');
    title(cond_lab{cond_ix});
    set(gca,'FontSize',14);
    
    % Misses
    subplot(4,numel(cond_lab),9+cond_ix); hold on;
    if strcmp(cond_lab{cond_ix},'Tar')
        scatter(1:numel(SBJs),n_miss_hit(:,cond_ix),sz,SBJ_colors,'filled');
        ylabel('# Missed Hit');
    else
        scatter(1:numel(SBJs),n_false_hit(:,cond_ix),sz,SBJ_colors,'filled');
        ylabel('# False Alarms');
    end
    line(xlim,[2 2],'Color','r');
    set(gca,'XLim',[0 numel(SBJs)+1]);
    set(gca,'XTick',1:numel(SBJs));
    set(gca,'XTickLabel',SBJs);
    xtickangle(45);
    xlabel('SBJ');
    title(cond_lab{cond_ix});
    set(gca,'FontSize',14);
end
