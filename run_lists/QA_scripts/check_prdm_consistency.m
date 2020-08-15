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
SBJs = {'EP07','EP08','EP09','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
    'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG09','EEG10','EEG12',...
    'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
    'EEG21','EEG22','EEG23','EEG24','EEG25','EEG26','EEG27','EEG28','EEG29','EEG30','EEG31'};
% 0 = good, 1 = suspect, 2 = likely toss, 3 = bad, 4 = low trial, 5 = oscillatory
status = [1 0 3 0 0 0 3 0 0 4 0 ...
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
%   EP06- removed for this analysis because of different prdm version
%   EP09
% SBJs = {'EP06','EP07','EP08','EP10','EP11','EP14','EP15','EP16','EP17','EP18','EP19',...
%            'EEG01','EEG02','EEG03','EEG04','EEG05','EEG06','EEG07','EEG08','EEG10','EEG12'};
% SBJs = {'EEG13','EEG14','EEG15','EEG16','EEG17','EEG18','EEG19','EEG20',...
%     'EEG21','EEG22','EEG23'};

%% Plot discrepancies
main_prdm = load([root_dir 'PRJ_Error_eeg/data/' SBJs{1} '/03_events/' SBJs{1} '_prdm_vars.mat']);
field_names = fieldnames(main_prdm);
for s = 2:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    
    prdm = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' SBJs{s} '_prdm_vars.mat']);
    fnames = fieldnames(prdm);
    if numel(fnames)~=numel(field_names)
        fprintf(2,'\t%s: # fields different!\n',SBJs{s});
    end
    if ~all(strcmp(fnames,field_names))
        fprintf(2,'\t%s: field names different!\n',SBJs{s});
    end
    for f = 1:numel(fnames)
        if any(strcmp(fnames{f},{'prdm_name','prdm_version'}))
            if ~strcmp(prdm.(fnames{f}),main_prdm.(fnames{f}))
                fprintf(2,'\t%s %s: %s vs. %s\n',SBJs{s},fnames{f},prdm.(fnames{f}),main_prdm.(fnames{f}));
            end
        elseif strcmp(fnames{f},'ITIs')
            tmp = 0;
        else
            if ~all(prdm.(fnames{f})==main_prdm.(fnames{f}))
                fprintf(2,'\t%s %s different!\n',SBJs{s},fnames{f});
            end
        end
    end
    
    clear SBJ_vars prdm
end

%% Plot individual fields: tol_lim, prdm_version, ITIs
for s = 1:numel(SBJs)
    prdm = load([root_dir 'PRJ_Error_eeg/data/' SBJs{s} '/03_events/' SBJs{s} '_prdm_vars.mat']);
    fprintf('%s = ',SBJs{s});
    disp(prdm.tol_lim);
end