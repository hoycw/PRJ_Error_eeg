%% Load BHV
SBJ = 'EEG14';
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

[bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJ '_behav.csv']);
%load([SBJ_vars.dirs.events SBJ '_behav_' proc_id '_final.mat']);

%% Plot SBJ Behavior Timeline
figure; hold on;
% RTs
plot(bhv.rt)
% RT bounds
line(xlim,[0.6 0.6],'Color','r');
line(xlim,[1.4 1.4],'Color','r');
% Tolerance
plot(bhv.tol+1,'Color','k')
plot(1-bhv.tol,'Color','k')

% Bad Feedback Losses
bf_l_idx = strcmp(bhv.fb,'L') & bhv.bad_fb;
scatter(find(bf_l_idx),bhv.rt(bf_l_idx),20,'r');
% Losses
gf_l_idx = strcmp(bhv.fb,'L') & ~bf_l_idx;
scatter(find(gf_l_idx),bhv.rt(gf_l_idx),40,'r','filled');

% Bad Feedback Wins
bf_w_idx = strcmp(bhv.fb,'W') & bhv.bad_fb;
scatter(find(bf_w_idx),bhv.rt(bf_w_idx),20,'g');
% Wins
gf_w_idx = strcmp(bhv.fb,'W') & ~bf_w_idx;
scatter(find(gf_w_idx),bhv.rt(gf_w_idx),40,'g','filled');

% Surprising
scatter(find(strcmp(bhv.fb,'S')),bhv.rt(strcmp(bhv.fb,'S')),75,'b','d');
