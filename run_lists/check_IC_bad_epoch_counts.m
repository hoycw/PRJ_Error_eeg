root_dir = '/Volumes/hoycw_clust/';
SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
    'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
    'EEG21','EEG22','EEG23'};
% SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
%            'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12'};
% SBJs = {'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
%     'EEG21','EEG22','EEG23'};
n_ic = zeros(size(SBJs));
n_be = zeros(size(SBJs));
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);

    training_ix = find(bhv.blk==0);
rt_low_ix   = find(bhv.rt <= proc.rt_bounds(1));
rt_high_ix  = find(bhv.rt >= proc.rt_bounds(2));
exclude_trials = unique(vertcat(bad_raw_trials, training_ix, rt_low_ix, rt_high_ix));
fprintf(2,'\tWarning: Removing %i trials (%i bad_raw, %i training, %i rts)\n', numel(exclude_trials),...
    numel(bad_raw_trials), numel(training_ix), numel(rt_low_ix)+numel(rt_high_ix));

    n_ic(s) = numel(SBJ_vars.ica_reject);
    n_be(s) = size(bad_epochs,1);
    clear SBJ_vars bad_epochs
end

%% 
scatter(n_ic,n_be);
xlabel('# ICs Rejected (n / 64)');
ylabel('# Raw Bad Epochs');