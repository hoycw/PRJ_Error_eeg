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
SBJs = {'EP06','EP07','EP08','EP09','EP10',...
    'EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
    'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
    'EEG21','EEG22','EEG23','EEG24','EEG25','EEG26','EEG27','EEG28','EEG29','EEG30','EEG31'};
% 0 = good, 1 = suspect, 2 = likely toss, 3 = bad, 4 = low trial, 5 = oscillatory
status = [3 1 0 3 0 ...
    0 0 3 0 0 4 0 ...
    0 3 0 0 0 0 4 0 4 0 1 ...
    0 1 0 0 0 1 0 0 ...
    0 0 0 0 0 0 0 1 0 0 2];
% status2 = zeros(size(SBJs));
% status2 = struct(...
%     'EP08','EP10','EP11','EP15','EP16','EP17','EP19',... % 0 good
%     'EP07',... % 1 suspect
%     'EP06','EP09','EP14',... % 3 bad
%     'EP18',... % 4 low trial
%     ...
%     'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12',...
%     'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
%     'EEG21','EEG22','EEG23','EEG24','EEG25','EEG26','EEG27','EEG28','EEG29','EEG30','EEG31'
% for s = 1:numel(SBJs)
%     status2.(SBJs{s}) = 
% end
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
n_ic_rej = zeros(size(SBJs));
n_bad_epochs = zeros(size(SBJs));
n_bad_epochs_tp = zeros(size(SBJs));
n_rt = zeros(size(SBJs));
n_bt = zeros(size(SBJs));
n_bc = zeros(size(SBJs));
n_bad_fb = zeros(size(SBJs));
perc_bad_fb = n_bad_fb;
n_good = zeros(size(SBJs));
quest_diff = zeros([numel(SBJs) numel(ideal_ans)]);
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    quest_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/questionare_vars/' SBJs{s} '_questionnaire_vars.m'];
    try
        eval(quest_vars_cmd);
        for q = 1:numel(ideal_ans)
            quest_diff(s,q) = abs(ideal_ans(q)-question_answers(q));
        end
    catch
    end
    
    [bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJs{s} '_behav.csv']);
    
    bad_epochs = fn_combine_raw_bad_epochs(SBJs{s});
    n_bad_points = 0;
    for b = 1:size(bad_epochs,1)
        n_bad_points = n_bad_points + numel(bad_epochs(b,1):bad_epochs(b,2));
    end
    
    n_bad_fb(s) = sum(bhv.bad_fb);
    perc_bad_fb(s) = sum(bhv.bad_fb)/numel(bhv.bad_fb);
    
    training_ix = find(bhv.blk==0);
    rt_low_ix   = find(bhv.rt <= proc.rt_bounds(1));
    rt_high_ix  = find(bhv.rt >= proc.rt_bounds(2));
%     exclude_trials = unique(vertcat(bad_raw_trials, training_ix, rt_low_ix, rt_high_ix));
%     fprintf(2,'\tWarning: Removing %i trials (%i bad_raw, %i training, %i rts)\n', numel(exclude_trials),...
%         numel(bad_raw_trials), numel(training_ix), numel(rt_low_ix)+numel(rt_high_ix));
    
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    n_good = numel(bhv.rt);
    
    n_ic_rej(s) = numel(SBJ_vars.ica_reject);
    n_bad_epochs(s) = size(bad_epochs,1);
    n_bad_epochs_tp(s) = n_bad_points;
    n_rt(s) = numel(unique(vertcat(rt_low_ix,rt_high_ix)));
    n_bt(s) = numel(SBJ_vars.trial_reject_ix);
    n_bc(s) = numel(SBJ_vars.ch_lab.bad);
    clear SBJ_vars bad_epochs bhv
end

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
scatter(1:numel(SBJs),n_bad_epochs,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# Bad Epochs');
set(gca,'FontSize',14);

subplot(7,1,3);
scatter(1:numel(SBJs),n_bad_epochs_tp,sz,SBJ_colors,'filled');
set(gca,'XTick',1:numel(SBJs));
set(gca,'XTickLabel',SBJs);
xtickangle(45);
xlabel('SBJ');
title('# Bad Epoch Time Points');
set(gca,'FontSize',14);

subplot(7,1,4);
scatter(1:numel(SBJs),n_ic_rej,sz,SBJ_colors,'filled');
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

% subplot(7,1,7);
% scatter(1:numel(SBJs),sum(n_bad_fb,2),sz,SBJ_colors,'filled');
% set(gca,'XTick',1:numel(SBJs));
% set(gca,'XTickLabel',SBJs);
% xtickangle(45);
% xlabel('SBJ');
% title('# bad feedback');
% set(gca,'FontSize',14);

%% 
scatter(n_ic_rej,n_bad_epochs);
xlabel('# ICs Rejected (n / 64)');
ylabel('# Raw Bad Epochs');