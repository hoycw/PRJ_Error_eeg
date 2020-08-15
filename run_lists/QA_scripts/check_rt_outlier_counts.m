root_dir = '/Volumes/hoycw_clust/';
proc_id  = 'eeg_full_ft';
SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
    'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
    'EEG21','EEG22','EEG23'};
proc_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/proc_vars/' proc_id '_vars.m'];
eval(proc_vars_cmd);
% SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
%            'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12'};
% SBJs = {'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
%     'EEG21','EEG22','EEG23'};
n_low  = zeros(size(SBJs));
n_high = zeros(size(SBJs));
acc_e  = zeros(size(SBJs));
acc_h  = zeros(size(SBJs));
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    [bhv] = fn_load_behav_csv([SBJ_vars.dirs.events SBJs{s} '_behav.csv']);
    
    e_idx = strcmp(bhv.cond,'easy');
    acc_e(s) = sum(bhv.hit(e_idx))/sum(e_idx);
    acc_h(s) = sum(bhv.hit(~e_idx))/sum(~e_idx);
    
    rt_low_ix   = find(bhv.rt <= proc.rt_bounds(1));
    rt_high_ix  = find(bhv.rt >= proc.rt_bounds(2));
    n_low(s)  = numel(rt_low_ix);
    n_high(s) = numel(rt_high_ix);
    fprintf(2,'%s: %i outliers (%i low; %i high); acc_e = %.03f, acc_h = %.03f\n', SBJs{s},...
        n_low(s)+n_high(s), n_low(s), n_high(s), acc_e(s), acc_h(s));
    
    clear SBJ_vars bhv
end
n_outliers = n_low+n_high;

%% 
figure;
subplot(2,1,1);
scatter(n_outliers,acc_e);
title('easy');
xlabel('# RT Outliers');
ylabel('Accuracy');

subplot(2,1,2);
scatter(n_outliers,acc_h);
title('hard');
xlabel('# RT Outliers');
ylabel('Accuracy');

%%
figure;
scatter(1-acc_e,acc_h);
xlabel('Easy Accuracy');
ylabel('Hard Accuracy');

