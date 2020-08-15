function QA_visualize
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

proc_id = 'eeg_full_ft';
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

%% Load data
% ICs tossed, bad trials tossed, channels tossed, bad epochs, bad RTs
%ideal_ans = [Easy Win, Easy Loss, Easy Surprise, Hard Win, Hard Loss, Hard
%Surprise]
ideal_ans = [7, 1, 5, 9, 3, 5];

%root_dir = '/Volumes/hoycw_clust/';
SBJs = {'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07', 'EEG08','EEG09','EEG10','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
    'EEG21','EEG22','EEG23','EEG24','EEG25','EEG26','EEG27','EEG28','EEG29','EEG30','EEG31'};
% 0 = good, 1 = suspect, 2 = likely toss, 3 = bad, 4 = low trial, 5 = oscillatory
status = [3 1 0 3 0 0 0 3 0 0 4 0 ...
    0 3 0 0 0 0 4 0 4 0 1 ...
    0 1 0 0 0 1 0 0 ...
    0 0 0 0 0 0 0 1 0 0 2];
colors = [...
    0 0 0; ...good- black
    1 0.5 0; ...suspect- orange
    1 0 1; ...likely toss - magenta
    1 0 0 ; ... bad- red
    0 1 0; ...low trial- green
    0 0 1]; % oscillatory- blue
SBJ_colors = zeros([numel(SBJs) 3]);
for s = 1:numel(SBJs)
    SBJ_colors(s,:) = colors(status(s)+1,:);
end
% Bad SBJ:
%   EP06
%   EP09
% SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
%            'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12'};
% SBJs = {'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
%     'EEG21','EEG22','EEG23'};

%%
qa_answers_all = zeros([numel(SBJs) numel(ideal_ans)]);
norm_crosscorr = zeros(size(SBJs));
outcome_diff = zeros(length(SBJs),3); %Win, Loss, Surprise
deviation_expected = zeros(length(SBJs),6);
for s = 1:numel(SBJs)
    display(SBJs{s})
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    quest_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/questionnaire_vars/' SBJs{s} '_questionnaire_vars.m'];
    try
        eval(quest_vars_cmd);
        for q = 1:numel(ideal_ans)
            deviation_expected(s,q) = abs(ideal_ans(q)-question_answers(q));
            qa_answers_all(s,q) = question_answers(q);
        end
    catch
    end
    for s = 1:numel(SBJs)
        temp = normxcorr2(ideal_ans,qa_answers_all(s,:));
        norm_crosscorr(s) = temp(length(ideal_ans));
        outcome_diff(s,2) = qa_answers_all(s,5) - qa_answers_all(s,2); %Hard Loss - Easy Loss
        outcome_diff(s,3) = qa_answers_all(s,6) - qa_answers_all(s,3); % Hard Surprise - Easy Surprise
        outcome_diff(s,1) = qa_answers_all(s,4) - qa_answers_all(s,1); % Hard Win - Easy Win
    end
 
end

%% 
sz = 75;
figure;
subplot(4,1,1);
scatter(1:numel(SBJs),norm_crosscorr,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(30);
%xlabel('SBJ');
title('Normalized Cross Correlation With Ideal');
set(gca,'FontSize',10);

subplot(4,1,2);
scatter(1:numel(SBJs), outcome_diff(:,1),sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(30);
%xlabel('SBJ');
yline(0)
title('Differences Between Wins In Hard - Easy Conditions');
set(gca,'FontSize',10);

subplot(4,1,3);
scatter(1:numel(SBJs), outcome_diff(:,2) ,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(30);
yline(0)
%xlabel('SBJ');
title('Differences Between Losses In Hard - Easy Conditions');
set(gca,'FontSize',10);

subplot(4,1,4);
scatter(1:numel(SBJs), outcome_diff(:,3),sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(30);
yline(0)
%xlabel('SBJ');
title('Differences Between Surprise Outcomes In Hard - Easy Conditions');
set(gca,'FontSize',10);
%WARNING: WILL ERROR FOR COLIN -- wasn't sure where in PRJ_Error to store
%them!
savefig('/Users/sheilasteiner/Desktop/scatter_plot_QA')
a = figure
violinplot_mod(deviation_expected, {'Easy Win', 'Easy Loss', 'Easy Surprise', 'Hard Win', 'Hard Loss', 'Hard Surprise'}, SBJs, status, SBJ_colors, 'Width', 0.45, 'ViolinAlpha', 0.2)
savefig(a, '/Users/sheilasteiner/Desktop/violin_plot_QA')
end
