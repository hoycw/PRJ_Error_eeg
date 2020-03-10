%% Load BHV
SBJ = 'EEG04';
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
bf_l_idx = bhv.hit==0 & bhv.bad_fb;
scatter(find(bf_l_idx),bhv.rt(bf_l_idx),20,'r');
% Losses
gf_l_idx = bhv.hit==0 & ~bf_l_idx;
scatter(find(gf_l_idx),bhv.rt(gf_l_idx),40,'r','filled');

