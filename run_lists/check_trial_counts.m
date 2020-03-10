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
[cond_lab,~,~,~,~] = fn_condition_label_styles('DifFB');
[eh_lab,~,~,~,~]   = fn_condition_label_styles('Dif');
n_trials_fb = zeros([numel(SBJs) numel(cond_lab)]);
n_trials_eh = zeros([numel(SBJs) numel(cond_lab)]);
for s = 1:numel(SBJs)
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Error_eeg/scripts/SBJ_vars/' SBJs{s} '_vars.m'];
    eval(SBJ_vars_cmd);
    load([SBJ_vars.dirs.events SBJs{s} '_behav_' proc_id '_final.mat']);
    
    fnames = fieldnames(bhv);
    for f = 1:numel(fnames)
        if ~strcmp(fnames{f},'blk')
            bhv.(fnames{f}) = bhv.(fnames{f})(bhv.blk~=1);
        end
    end
    bhv.blk = bhv.blk(bhv.blk~=1);
    
    % Dif FB
    cond_idx = fn_condition_index(cond_lab,bhv);
    for cond_ix = 1:numel(cond_lab)
        n_trials_fb(s,cond_ix) = sum(cond_idx==cond_ix);
    end
    
    % Dif
    cond_idx = fn_condition_index(eh_lab,bhv);
    for cond_ix = 1:numel(eh_lab)
        n_trials_eh(s,cond_ix) = sum(cond_idx==cond_ix);
    end
    
    clear SBJ_vars bhv cond_idx
end

%% 
sz = 75;
figure;
for cond_ix = 1:numel(cond_lab)
    subplot(numel(cond_lab),1,cond_ix); hold on;
    scatter(1:numel(SBJs),n_trials_fb(:,cond_ix),sz,SBJ_colors,'filled');
    line(xlim,[30 30],'Color','r');
    ylims = ylim;
    set(gca,'YLim',[0 ylims(2)]);
    set(gca,'XTick',1:numel(SBJs));
    set(gca,'XTickLabel',SBJs);
    xtickangle(45);
    xlabel('SBJ');
    ylabel('# Trials');
    title(cond_lab{cond_ix});
    set(gca,'FontSize',14);
end

figure;
for cond_ix = 1:numel(eh_lab)
    subplot(numel(eh_lab),1,cond_ix); hold on;
    scatter(1:numel(SBJs),n_trials_eh(:,cond_ix),sz,SBJ_colors,'filled');
%     line(xlim,[30 30],'Color','r');
    set(gca,'XTick',1:numel(SBJs));
    set(gca,'XTickLabel',SBJs);
    xtickangle(45);
    xlabel('SBJ');
    ylabel('# Trials');
    title(eh_lab{cond_ix});
    set(gca,'FontSize',14);
end
