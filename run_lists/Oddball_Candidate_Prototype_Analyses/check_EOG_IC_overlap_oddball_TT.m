% Check overlap in EOG ICs from Oddball trials and Target Time trials
eogs_missing = false(size(SBJs));
eog_diff     = false(size(SBJs));
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/02_preproc/' SBJs{s} '_odd_full_ft_02a.mat'],'heog_ics','veog_ics');
    tmp_h = heog_ics;
    tmp_v = veog_ics;
    load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/02_preproc/' SBJs{s} '_eeg_full_ft_02a.mat'],'heog_ics','veog_ics');
    all_bad = unique([tmp_h tmp_v heog_ics veog_ics]);
    reject_match = false(size(all_bad));
    for b = 1:numel(all_bad)
        reject_match(b) = any(all_bad(b)==SBJ_vars.ica_reject);
    end
    if all(reject_match)
        fprintf('%s: all EOG ICs in reject_ix\n',SBJs{s});
    else
        fprintf(2,'%s: EOG ICs missing from reject_ix~\n',SBJs{s});
        eogs_missing(s) = true;
    end
    if ~numel(tmp_h)==numel(heog_ics) || ~all(tmp_h==heog_ics)
        fprintf(2,'%s: H not same; odd = %s; eeg = %s\n',SBJs{s},num2str(tmp_h),num2str(heog_ics));
        eog_diff(s) = true;
    end
    if ~numel(tmp_v)==numel(veog_ics) || ~all(tmp_v==veog_ics)
        fprintf(2,'%s: V not same; odd = %s; eeg = %s\n',SBJs{s},num2str(tmp_v),num2str(veog_ics));
        eog_diff(s) = true;
    end
end

fprintf('%d / %d SBJs missing EOG ICs from SBJ_vars.ica_reject\n',sum(eogs_missing),numel(SBJs));
fprintf('SBJs with missing: '); disp(SBJs(eogs_missing)');
fprintf('%d / %d SBJs mismatch between oddball and TT EOG ICs\n',sum(eog_diff),numel(SBJs));
fprintf('SBJs with mismatch odd/TT EOG ICs: '); disp(SBJs(eog_diff)');