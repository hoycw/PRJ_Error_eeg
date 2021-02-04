%% Plot Auditory Salience Metrics for OB and TT stimuli
if exist('/home/knight/','dir');root_dir='/home/knight/';app_dir=[root_dir 'PRJ_Error_eeg/Apps/'];
elseif exist('/Users/sheilasteiner/','dir'); root_dir='/Users/sheilasteiner/Desktop/Knight_Lab/';app_dir='/Users/sheilasteiner/Downloads/fieldtrip-master/';
else root_dir='/Volumes/hoycw_clust/';app_dir='/Users/colinhoy/Code/Apps/';end

ob_color   = [166 86 40]./255;      % brown
ob_style   = '--';
tt_color   = [55 126 184]./255;     % blue
tt_style   = '-';
win_color  = [77 175 74]./255;      % green
loss_color = [228 26 28]./255;      % red
n_bins     = 25;
face_alpha = 0.6;
fig_ftype  = 'png';

%% Load data and set up stats
[reg_lab, reg_names, ~, ~, ~] = fn_regressor_label_styles('AudSal');

% Load auditory salience features
aud_sal = fn_load_auditory_salience([root_dir 'PRJ_Error_eeg/scripts/auditory_salience_output.csv']);
metrics = fieldnames(aud_sal);
metrics = metrics(~strcmp(metrics,'sound'));

% Sort by feedback and paradigm
win_ix  = find(strcmp(aud_sal.sound,'new_win_sound.wav'));
loss_ix = find(strcmp(aud_sal.sound,'new_loss_sound.wav'));
surp_ix = find(~cellfun(@isempty,strfind(aud_sal.sound,'.wav')));
surp_ix = surp_ix(surp_ix~=win_ix & surp_ix~=loss_ix);
ob_ix   = find(~cellfun(@isempty,strfind(aud_sal.sound,'.WAV')));

% Get range of data
ranges = nan([numel(metrics) 2]);
ob_means = nan(size(metrics));
ob_vars  = nan(size(metrics));
tt_means = nan(size(metrics));
tt_vars  = nan(size(metrics));
for m = 1:numel(metrics)
    ranges(m,:) = [min(aud_sal.(metrics{m})) max(aud_sal.(metrics{m}))];
    ob_means(m) = mean(aud_sal.(metrics{m})(ob_ix));
    ob_vars(m)  = std(aud_sal.(metrics{m})(ob_ix));
    tt_means(m) = mean(aud_sal.(metrics{m})(surp_ix));
    tt_vars(m)  = std(aud_sal.(metrics{m})(surp_ix));
end

%% Plot Summary Statistics
% Create figure directory
fig_dir = [root_dir 'PRJ_Error_eeg/results/model_predictions/auditory_salience/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Plot design matrix
fig_name = 'auditory_salience_properties';
figure('Name',fig_name,'units','normalized','outerposition',[0 0 1 1]);

for m = 1:numel(metrics)
    subplot(2,3,m); hold on;
    % Plot Histograms
    ob_hist = histogram(aud_sal.(metrics{m})(ob_ix),linspace(ranges(m,1),ranges(m,2),n_bins),...
        'FaceColor',ob_color,'FaceAlpha',face_alpha);
    tt_hist = histogram(aud_sal.(metrics{m})(surp_ix),linspace(ranges(m,1),ranges(m,2),n_bins),...
        'FaceColor',tt_color,'FaceAlpha',face_alpha);
    
    % Plot Means
    ob_line = line([ob_means(m) ob_means(m)],ylim,'Color','k','LineWidth',3,'LineStyle',ob_style);
    tt_line = line([tt_means(m) tt_means(m)],ylim,'Color','k','LineWidth',3,'LineStyle',tt_style);
    
    % Plot Win and Loss
    win_line = line([aud_sal.(metrics{m})(win_ix) aud_sal.(metrics{m})(win_ix)],...
        ylim,'Color',win_color,'LineWidth',2);
    loss_line = line([aud_sal.(metrics{m})(loss_ix) aud_sal.(metrics{m})(loss_ix)],...
        ylim,'Color',loss_color,'LineWidth',2);
    
    % Properties
    xlim(ranges(m,:));
    legend([ob_hist, tt_hist, ob_line, tt_line, win_line, loss_line],...
        {['Oddball (n=' num2str(numel(ob_ix)) ')'],['Target Time (n=' num2str(numel(surp_ix)) ')'],...
        ['mean(OB) = ' num2str(ob_means(m),'%.2f') '+/-' num2str(ob_vars(m),'%.2f')],...
        ['mean(TT) = ' num2str(tt_means(m),'%.2f') '+/-' num2str(tt_vars(m),'%.2f')],...
        ['Win Sound=' num2str(aud_sal.(metrics{m})(win_ix),'%.2f')],...
        ['Loss Sound=' num2str(aud_sal.(metrics{m})(loss_ix),'%.2f')]},'Location','northeast');
    title(reg_names{strcmp(reg_lab,metrics{m})});
    set(gca,'FontSize',16);
end

%% Save Figure
fig_fname = [fig_dir fig_name '.' fig_ftype];
fprintf('Saving %s\n',fig_fname);
% Ensure vector graphics if saving
if any(strcmp(fig_ftype,{'svg','eps'}))
    set(gcf, 'Renderer', 'painters');
end
saveas(gcf,fig_fname);
