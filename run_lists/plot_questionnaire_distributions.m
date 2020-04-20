if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

proc_id = 'eeg_full_ft';
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

%% Question Parameters
q_lab = {'Easy Win','Easy Loss','Easy Surp','Hard Win','Hard Loss','Hard Surp'};
ideal_ans = [7, 1, 5, 9, 3, 5];
bounds = [6 nan; nan 4; 4 6; 6 nan; nan 4; 4 6];

%% Load data
SBJ_id   = 'EEGall';
good_id  = 'goodEEG';
maybe_id = 'maybe';
bad_id   = 'bad';

root_dir = '/Volumes/hoycw_clust/';
SBJs = fn_load_SBJ_list(SBJ_id);
good_SBJs = fn_load_SBJ_list(good_id);
maybe_SBJs = fn_load_SBJ_list(maybe_id);
bad_SBJs = fn_load_SBJ_list(bad_id);

% SBJs = {'EP06','EP07','EP08','EP09','EP10',...
%     'EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
%     'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12',...
%     'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
%     'EEG21','EEG22','EEG23','EEG24','EEG25','EEG26','EEG27','EEG28','EEG29','EEG30','EEG31'};
% status = [3 1 0 3 0 ...
%     0 0 3 0 0 4 0 ...
%     0 3 0 0 0 0 4 0 4 0 1 ...
%     0 1 0 0 0 1 0 0 ...
%     0 0 0 0 0 0 0 1 0 0 2];

%% Assign labels/colors
% 0 = good, 1 = suspect, 2 = likely toss, 3 = bad, 4 = low trial, 5 = oscillatory
status = zeros(size(SBJs));
for s = 1:numel(SBJs)
    if any(strcmp(SBJs{s},maybe_SBJs))
        status(s) = 2;
    elseif any(strcmp(SBJs{s},bad_SBJs))
        status(s) = 3;
    elseif any(strcmp(SBJs{s},good_SBJs))
        status(s) = 1;
    else
        error([SBJs{s} ' not found!']);
    end
end
colors = [...
    0 0 0; ...good- black
    1 0.5 0; ...suspect- orange
    1 0 0]; % bad- red
SBJ_colors = zeros([numel(SBJs) 3]);
for s = 1:numel(SBJs)
    SBJ_colors(s,:) = colors(status(s),:);
end
% Bad SBJ:
%   EP06
%   EP09
% SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
%            'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12'};
% SBJs = {'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
%     'EEG21','EEG22','EEG23'};

%%
answers = zeros([numel(SBJs) numel(ideal_ans)]);
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    quest_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/questionnaire_vars/' SBJs{s} '_questionnaire_vars.m'];
    try
        eval(quest_vars_cmd);
        for q = 1:numel(q_lab)
            answers(s,q) = question_answers(q);
        end
    catch
        fprintf(2,'%s no answers found!\n',SBJs{s});
    end
    
    clear SBJ_vars quest_vars
end

%% Plot Answers
bins = 1:10;
sz = 75;
figure;
for q = 1:numel(q_lab)
    subplot(1,numel(q_lab),q); hold on;
    scatter(answers(:,q),1:numel(SBJs),sz,SBJ_colors,'filled');
    line([ideal_ans(q) ideal_ans(q)],ylim,'Color','k','LineWidth',2);
    line([bounds(q,1) bounds(q,1)],ylim,'Color','r');
    line([bounds(q,2) bounds(q,2)],ylim,'Color','r');
    xlim([0 10]);
    set(gca,'YTick',1:numel(SBJs));
    set(gca,'YTickLabel',SBJs);
    ytickangle(45);
    ylabel('SBJ');
    xlabel('bad <- rating -> good')
    title(q_lab{q});
    set(gca,'FontSize',14);
end
