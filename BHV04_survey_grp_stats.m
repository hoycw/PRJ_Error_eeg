function BHV04_survey_grp_stats(SBJ_id, plot_fig, save_fig, fig_ftype)
%% Stats and plotting of group level survey responses
% INPUTS:
%   SBJ_id [str] - name of list of SBJs
%   plot_fig [0/1] - binary flag to plot data
%   save_fig [0/1] - binary flag to save figure
%   fig_ftype [str] - file extension to save fig (e.g., 'png','svg', etc.)
% OUTPUTS:
%   optionally saves survey data summary figure

%% Set up paths
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

addpath([root_dir 'PRJ_Error_eeg/scripts/']);
addpath([root_dir 'PRJ_Error_eeg/scripts/utils/']);

%% Processing pipeline
proc_id = 'eeg_full_ft';
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

SBJs = fn_load_SBJ_list(SBJ_id);
[cond_lab, cond_names, cond_colors, ~, cond_mrkrs] = fn_condition_label_styles('DifFB');

% Define ideal answers and cut offs (for visualization only)
%   "How would you feel about this feedback?"
%   1 = "Terrible", 5 = "I don't care", 9 = "Great"
ideal_ans = [7, 5, 1, 9, 5, 3];
bounds = [6 nan; 4 6; nan 4; 6 nan; 4 6; nan 4];

%% Load data
answers = nan([numel(SBJs) numel(cond_lab)]);
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    quest_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/questionnaire_vars/' SBJs{s} '_questionnaire_vars.m'];
    
    try
        eval(quest_vars_cmd);
        for cond_ix = 1:numel(cond_lab)
            % Subtract 5 to mean center data on indifference point
            answers(s,strcmp(cond_lab,question_conditions{cond_ix})) = question_answers(cond_ix)-5;
        end
    catch
        fprintf(2,'%s no answers found!\n',SBJs{s});
    end
    
    clear SBJ_vars quest_vars question_answers question_conditions free_response
end

%% Print descriptive statistics
%   test each answer for valence via t-test against zero ("I don't care" neutral)
pval = nan(size(cond_lab));
tval = nan(size(cond_lab));
for cond_ix = 1:numel(cond_lab)
    [~,pval(cond_ix),~,stats] = ttest(answers(:,cond_ix));
    tval(cond_ix) = stats.tstat;
    fprintf('%s [mean +/- SD] = %.3f +/- %.3f, t(%d) = %.5f, p = %.5f\n',cond_lab{cond_ix},...
        nanmean(answers(:,cond_ix)-5),nanstd(answers(:,cond_ix)),stats.df,tval(cond_ix),pval(cond_ix));
end

%% Plot Data
if plot_fig
    % Plot Group Answers
    fig_name = [SBJ_id '_questionnaire_hist'];
    figure('Name',fig_name);
    for cond_ix = 1:numel(cond_lab)
        subplot(numel(cond_lab),1,cond_ix); hold on;
        title([cond_names{cond_ix} ': p=' num2str(pval(cond_ix),'%.5f')]);
        histogram(answers(:,cond_ix),'FaceColor',cond_colors{cond_ix});
        line([nanmean(answers(:,cond_ix)) nanmean(answers(:,cond_ix))], ylim,...
            'Color', 'k', 'LineWidth', 2);
        xlim([-6 6]);
    end
    
    if save_fig
        fig_dir = [root_dir 'PRJ_Error_eeg/results/questionnaire/'];
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end

    % Plot Single SBJ Answers
    bins = 1:10;
    sz = 75;
    fig_name = [SBJ_id '_questionnaire_scatter'];
    figure('Name',fig_name);
    for cond_ix = 1:numel(cond_lab)
        subplot(1,numel(cond_lab),cond_ix); hold on;
        scatter(answers(:,cond_ix),1:numel(SBJs),sz,'filled');
        line([ideal_ans(cond_ix) ideal_ans(cond_ix)]-5,ylim,'Color','k','LineWidth',2);
        line([bounds(cond_ix,1) bounds(cond_ix,1)]-5,ylim,'Color','r');
        line([bounds(cond_ix,2) bounds(cond_ix,2)]-5,ylim,'Color','r');
        xlim([-6 6]);
        set(gca,'YTick',1:numel(SBJs));
        set(gca,'YTickLabel',SBJs);
        ytickangle(45);
        ylabel('SBJ');
        xlabel('bad <- rating -> good')
        title(cond_names{cond_ix});
        set(gca,'FontSize',14);
    end
    
    if save_fig
        fig_dir = [root_dir 'PRJ_Error_eeg/results/questionnaire/'];
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
end

end