if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

proc_id = 'eeg_full_ft';
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);

%% Load data
% ICs tossed, bad trials tossed, channels tossed, bad epochs, bad RTs

ideal_ans = [7, 1, 5, 9, 3, 5];

root_dir = '/Volumes/hoycw_clust/';
SBJs = {'EP06','EP07','EP08','EP09','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
    'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
    'EEG21','EEG22','EEG23','EEG24','EEG25','EEG26','EEG27','EEG28','EEG29','EEG30','EEG31'};
% 0 = good, 1 = suspect, 2 = likely toss, 3 = bad, 4 = low trial, 5 = oscillatory
status = [2 1 0 2 0 0 0 3 5 0 4 0 ...
    0 3 0 0 0 0 4 0 4 0 4 ...
    0 0 0 5 0 1 0 0 ...
    5 0 0 0 0 0 0 1 5 0 4];
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
bhvs = cell(size(SBJs));
quest_SBJs = false(size(SBJs));
sbj_ans = zeros([numel(SBJs) numel(ideal_ans)]);
quest_diff = zeros([numel(SBJs) numel(ideal_ans)]);
for s = 1:numel(SBJs)
    fprintf('\tLoading %s\n',SBJs{s});
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    quest_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/questionare_vars/' SBJs{s} '_questionnaire_vars.m'];
    try
        eval(quest_vars_cmd);
        for q = 1:numel(ideal_ans)
            quest_diff(s,q) = abs(ideal_ans(q)-question_answers{q});
        end
        quest_SBJs(s) = true;
    catch
    end
    
    bhvs{s} = load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat'],'bhv');
    
    clear SBJ_vars bad_epochs bhv
end
qSBJs = SBJs(quest_SBJs);

%% 
sz = 75;
figure;
subplot(7,1,1);
scatter(1:numel(SBJs),n_bc,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# bad channels');
set(gca,'FontSize',14);

subplot(7,1,2);
scatter(1:numel(SBJs),n_be,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# Bad Epochs');
set(gca,'FontSize',14);

subplot(7,1,3);
scatter(1:numel(SBJs),n_be_tp,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# Bad Epoch Time Points');
set(gca,'FontSize',14);

subplot(7,1,4);
scatter(1:numel(SBJs),n_ic,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# ICs Tossed');
set(gca,'FontSize',14);

subplot(7,1,5);
scatter(1:numel(SBJs),n_rt,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# RTs Tossed');
set(gca,'FontSize',14);

subplot(7,1,6);
scatter(1:numel(SBJs),n_bt,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# Bad Trials (noise)');
set(gca,'FontSize',14);

subplot(7,1,7);
scatter(1:numel(SBJs),sum(quest_diff,2),sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('??? Difference from Ideal');
set(gca,'FontSize',14);
% 
% subplot(8,1,7);
% scatter(1:numel(SBJs),n_good,sz,SBJ_colors,'filled');
% set(gca,'XTick',1:numel(SBJs));
% set(gca,'XTickLabel',SBJs);
% xtickangle(45);
% xlabel('SBJ');
% title('# good RTs');
% set(gca,'FontSize',14);

%% 
scatter(n_ic,n_be);
xlabel('# ICs Rejected (n / 64)');
ylabel('# Raw Bad Epochs');