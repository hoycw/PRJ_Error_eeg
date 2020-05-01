%% Model Predictions
cond_lab = {'Likely Win', 'Unlikely Neu', 'Unlikely Loss', 'Unlikely Win', 'Unlikely Neu', 'Likely Loss'};
% cond_lab = {'Easy Win', 'Easy Neu', 'Easy Loss', 'Hard Win', 'Hard Neu', 'Hard Loss'};

% Categorical Reward
exp  = [1 1 1 -1 -1 -1];
out  = [1 0 -1 1 0 -1];
mag  = [1 0 1 1 0 1];
interaction = exp .* out;

% Continuous Reward Prediction Error
exp_rew = [0.6 0.6 0.6 -0.6 -0.6 -0.6];
rpe     = out - exp_rew;
urpe    = abs(rpe);

% Categorical Likelihood
out_cat  = [1 0 0 0 0 1];
% Continuous Likelihood
out_freq = [0.75 0.1 0.15 0.15 0.1 0.75];
out_rare = 1-out_freq;

%% Plot Models
fig_name  = 'model_prediction_comparison';
fig_ftype = 'png';
figure('Name',fig_name,'units','normalized','OuterPosition',[0 0 1 1]);

signed_color = [209 151 105]./255;
signed_style = '-';
signed_mark  = 's';

unsigned_color = [118 160 156]./255;
unsigned_style = ':';
unsigned_mark  = 'o';

lik_color = [152 78 163]./255;
lik_style = '--';
lik_mark  = '*';

alt_color = 'k';    % categorical salience/likelihood/frequency
odd_color = 'r';    % meaningless interaction
alt_style = '-';
alt_mark  = 's';

cond_idx = 1:6;
easy_idx = 1:3;
hard_idx = 4:6;

cond_idx_nn = [1 3 4 6];
easy_idx_nn = [1 3];
hard_idx_nn = [4 6];

offset     = 0.3;
tick_angle = -20;
width      = 2;

%% ======================= Without Neutral Trials =======================
%--------------------- Categorical ---------------------
subplot(2,3,1); hold on;
title('Categorical Outcomes (No Neu.)');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx_nn);
set(gca,'XTickLabels',cond_lab(cond_idx_nn));
xtickangle(tick_angle);
xlim([0 7]);

yyaxis left
ylabel('Reward Category');
tmp1 = plot(easy_idx_nn, out(easy_idx_nn), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
plot(hard_idx_nn, out(hard_idx_nn), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
ylim([-1.1 1.1]);
yticks([-1 1]);
yticklabels({'(Loss) -1','(Win) 1'});
set(gca,'YColor',signed_color);

yyaxis right
ylabel('Likelihood Category');
tmp2 = plot(easy_idx_nn+offset, interaction(easy_idx_nn), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
plot(hard_idx_nn+offset, interaction(hard_idx_nn), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
ylim([-1.1 1.1]);
yticks([-1 1]);
yticklabels({'1 (Unlikely)','-1 (Likely)'});
set(gca,'YColor',unsigned_color);
legend([tmp1 tmp2], {'Outcome','Interaction'}, 'Location','southwest');

%--------------------- Continuous RPE ---------------------
subplot(2,3,2); hold on;
title('Continuous Prediction Error (No Neu.)');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx_nn);
set(gca,'XTickLabels',cond_lab(cond_idx_nn));
xtickangle(tick_angle);
xlim([0 7]);

yyaxis left
ylabel('Reward Prediction Error');
tmp1 = plot(easy_idx_nn, rpe(easy_idx_nn), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
plot(hard_idx_nn, rpe(hard_idx_nn), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
ylim([-2.1 2.1]);
yticks([-2:1:2]);
set(gca,'YColor',signed_color);

yyaxis right
ylabel('PE Magnitude/Rareness');
tmp2 = plot(easy_idx_nn+offset, urpe(easy_idx_nn), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
plot(hard_idx_nn+offset, urpe(hard_idx_nn), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
ylim([-2.1 2.1]);
yticks([-2:1:2]);
set(gca,'YColor',unsigned_color);

tmp3 = plot(easy_idx_nn+2*offset, out_rare(easy_idx_nn), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
plot(hard_idx_nn+2*offset, out_rare(hard_idx_nn), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
legend([tmp1 tmp2 tmp3], {'RPE','RPE Magn.','RPE Rare'}, 'Location','southeast');

%--------------------- Salience ---------------------
subplot(2,3,3); hold on;
title('Pure Salience (No Neu.)');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx_nn);
set(gca,'XTickLabels',cond_lab(cond_idx_nn));
xtickangle(tick_angle);
xlim([0 7]);

yyaxis left
ylabel('Salience Category');
tmp1 = plot(easy_idx_nn, out_cat(easy_idx_nn), 'Color', alt_color, 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
plot(hard_idx_nn, out_cat(hard_idx_nn), 'Color', alt_color, 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'(Unlikely) 0','(Likely) 1'});
set(gca,'YColor',alt_color);

yyaxis right
ylabel('Proportion of Trials');
tmp2 = plot(easy_idx_nn+offset, out_freq(easy_idx_nn), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
plot(hard_idx_nn+offset, out_freq(hard_idx_nn), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
ylim([-0.1 1.1]);
set(gca,'YColor',unsigned_color);
legend([tmp1 tmp2], {'Categorical Likelihood','Continuous Likelihood'}, 'Location','north');

%% ======================= With Neutral Trials =======================
%--------------------- Categorical ---------------------
subplot(2,3,4); hold on;
title('Categorical Outcomes');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx);
set(gca,'XTickLabels',cond_lab);
xtickangle(tick_angle);
xlim([0 7]);

yyaxis left
ylabel('Reward Value');
tmp1 = plot(easy_idx, out(easy_idx), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
plot(hard_idx, out(hard_idx), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
ylim([-1.1 1.1]);
yticks([-1 0 1]);
yticklabels({'(Loss) -1','(Neu) 0','(Win) 1'});
set(gca,'YColor',signed_color);

yyaxis right
ylabel('Reward Magnitude');
tmp2 = plot(easy_idx+offset, mag(easy_idx), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
plot(hard_idx+offset, mag(hard_idx), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
ylim([-1.1 1.1]);
yticks([-1 0 1]);
yticklabels({'\color{red}-1 (???)','0 (Low)','1 (High)'})
set(gca,'YColor',unsigned_color);

tmp3 = plot(easy_idx+2*offset, interaction(easy_idx), 'Color', 'r', 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
plot(hard_idx+2*offset, interaction(hard_idx), 'Color', 'r', 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
legend([tmp1 tmp2 tmp3], {'Outcome','Magnitude','Interaction'}, 'Location','southwest');

%--------------------- Continuous RPE ---------------------
subplot(2,3,5); hold on;
title('Continuous Prediction Error');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx);
set(gca,'XTickLabels',cond_lab);
xtickangle(tick_angle);
xlim([0 7]);

yyaxis left
ylabel('Reward Prediction Error');
tmp1 = plot(easy_idx, rpe(easy_idx), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
plot(hard_idx, rpe(hard_idx), 'Color', signed_color, 'LineStyle', signed_style, 'LineWidth', width, 'Marker', signed_mark);
ylim([-2.1 2.1]);
yticks([-2:1:2]);
set(gca,'YColor',signed_color);

yyaxis right
ylabel('PE Magnitude/Rareness');
tmp2 = plot(easy_idx+offset, urpe(easy_idx), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
plot(hard_idx+offset, urpe(hard_idx), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
ylim([-2.1 2.1]);
yticks([-2:1:2]);
set(gca,'YColor',unsigned_color);

tmp3 = plot(easy_idx+2*offset, out_rare(easy_idx), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
plot(hard_idx+2*offset, out_rare(hard_idx), 'Color', lik_color, 'LineStyle', lik_style, 'LineWidth', width, 'Marker', lik_mark);
legend([tmp1 tmp2 tmp3], {'RPE','RPE Magn.','RPE Rare'}, 'Location','southeast');

%--------------------- Salience ---------------------
subplot(2,3,6); hold on;
title('Pure Salience');
set(gca,'FontSize',16);
set(gca,'XTick',cond_idx);
set(gca,'XTickLabels',cond_lab);
xtickangle(tick_angle);
xlim([0 7]);

yyaxis left
ylabel('Salience Category');
tmp1 = plot(easy_idx, out_cat(easy_idx), 'Color', alt_color, 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
plot(hard_idx, out_cat(hard_idx), 'Color', alt_color, 'LineStyle', alt_style, 'LineWidth', width, 'Marker', alt_mark);
ylim([-0.1 1.1]);
yticks([0 1]);
yticklabels({'(Unlikely) 0','(Likely) 1'});
set(gca,'YColor',alt_color);

yyaxis right
ylabel('Proportion of Trials');
tmp2 = plot(easy_idx+offset, out_freq(easy_idx), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
plot(hard_idx+offset, out_freq(hard_idx), 'Color', unsigned_color, 'LineStyle', unsigned_style, 'LineWidth', width, 'Marker', unsigned_mark);
ylim([-0.1 1.1]);
set(gca,'YColor',unsigned_color);
legend([tmp1 tmp2], {'Categorical Likelihood','Continuous Likelihood'}, 'Location','north');

%% Save figure
fig_dir = [root_dir 'PRJ_Error_eeg/results/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

fig_fname = [fig_dir fig_name '.' fig_ftype];
% Commented out because screen ratios aren't right automatically
% fprintf('Saving %s\n',fig_fname);
% saveas(gcf,fig_fname);
